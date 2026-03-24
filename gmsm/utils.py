
import cobra
import datetime
import logging
import math
import os
import pickle
import sys
import subprocess
from os.path import abspath, dirname, getmtime, isfile, join, split

_ELF_MAGIC = b"\x7fELF"
_MACHO_MAGICS = {
    b"\xfe\xed\xfa\xce",
    b"\xce\xfa\xed\xfe",
    b"\xfe\xed\xfa\xcf",
    b"\xcf\xfa\xed\xfe",
    b"\xca\xfe\xba\xbe",
    b"\xbe\xba\xfe\xca",
}

_TEMPLATE_SCORE_DEFAULTS = {
    'template_ani_weight': 0.7,
    'template_af_weight': 0.3,
    'template_diamond_hit_weight': 0.05,
    'template_diamond_identity_weight': 0.95,
    'template_bbh_template_weight': 0.5,
    'template_bbh_target_weight': 0.5,
    'template_coarse_weight': 0.95,
    'template_rerank_weight': 0.05,
}

_TEMPLATE_SCORE_PAIRS = (
    ('template_ani_weight', 'template_af_weight'),
    ('template_diamond_hit_weight', 'template_diamond_identity_weight'),
    ('template_bbh_template_weight', 'template_bbh_target_weight'),
    ('template_coarse_weight', 'template_rerank_weight'),
)


def setup_logging(run_ns):
    if run_ns.verbose:
        log_level = logging.INFO
    elif run_ns.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING

    logging.basicConfig(format='%(levelname)s: %(message)s', level=log_level)

    if run_ns.verbose or run_ns.debug:
        os.makedirs(run_ns.outputfolder, exist_ok=True)
        logger = logging.getLogger('')
        fomatter = logging.Formatter(
                '[%(levelname)s|%(filename)s:%(lineno)s] > %(message)s')
        fh = logging.FileHandler(
                os.path.join(run_ns.outputfolder, 'gmsm.log'), mode = 'w')
        fh.setFormatter(fomatter)
        logger.setLevel(log_level)
        logger.addHandler(fh)


def get_version():
    import gmsm
    version = gmsm.__version__

    return version


def get_git_log():
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        out, err, retcode = execute(args)
        if isinstance(out, bytes):
            return out.decode('utf-8', errors='replace').strip()
        return out.strip()
    except OSError:
        pass
    return""


def load_legacy_pickle(path):
    """Load old Python pickles shipped with GMSM across Python 3 versions."""
    from gmsm import runtime_assets

    resolved_path = runtime_assets.resolve_runtime_asset_path(path)
    if not os.path.isfile(resolved_path) or runtime_assets.is_lfs_pointer_file(resolved_path):
        raise FileNotFoundError(runtime_assets.runtime_asset_missing_message(path))

    with open(resolved_path, 'rb') as handle:
        raw = handle.read()

    try:
        return pickle.loads(raw)
    except (UnicodeDecodeError, pickle.UnpicklingError):
        normalized = raw.replace(b'\r\n', b'\n')
        return pickle.loads(normalized, encoding='latin1')


def load_legacy_cobra_pickle(path):
    """Load legacy COBRA pickles whose optlang solver state predates current fields."""
    try:
        import optlang.glpk_interface as glpk_interface
    except ImportError:
        return load_legacy_pickle(path)

    original_setstate = glpk_interface.Configuration.__setstate__

    def _patched_setstate(self, state):
        if "tolerances" not in state:
            state = dict(state)
            state["tolerances"] = {}
        return original_setstate(self, state)

    glpk_interface.Configuration.__setstate__ = _patched_setstate
    try:
        loaded = load_legacy_pickle(path)
    finally:
        glpk_interface.Configuration.__setstate__ = original_setstate

    if isinstance(loaded, cobra.Model):
        ensure_modern_cobra_attrs(loaded)
    return loaded


def ensure_modern_cobra_attrs(model):
    """Backfill attributes expected by current COBRApy on legacy pickled objects."""
    objects = [model]
    objects.extend(model.reactions)
    objects.extend(model.metabolites)
    objects.extend(model.genes)
    if hasattr(model, "groups"):
        objects.extend(model.groups)

    for obj in objects:
        if not hasattr(obj, "_annotation"):
            obj._annotation = {}


def check_input_options(run_ns):
    input_file = getattr(run_ns, 'input', None)
    ec_file = getattr(run_ns, 'ec_file', None)
    pmr_generation = getattr(run_ns, 'pmr_generation', False)
    smr_generation = getattr(run_ns, 'smr_generation', False)
    comp = getattr(run_ns, 'comp', None)
    template_topk = getattr(run_ns, 'template_topk', 3)
    template_genome_bank = getattr(run_ns, 'template_genome_bank', None)
    template_rerank_topn = getattr(run_ns, 'template_rerank_topn', 3)
    template_recommendation_only = getattr(run_ns, 'template_recommendation_only', False)

    if not input_file:
        logging.warning("Provide input file via ('-i')")
        sys.exit(1)

    if not ec_file and \
            not pmr_generation and \
            not smr_generation and \
            not comp:
                logging.warning("Select one of the options: '-p' or '-s'")
                sys.exit(1)

    if comp:
        if not pmr_generation:
            logging.warning(
                    "Primary metabolic modeling option ('-p') should also be selected")
            sys.exit(1)

    if ec_file:
        if not pmr_generation:
            logging.warning(
                    "Primary metabolic modeling option ('-p') should also be selected")
            sys.exit(1)

    if template_topk is not None and int(template_topk) < 1:
        logging.warning("Template recommendation requires '--template-topk' to be at least 1")
        sys.exit(1)

    if template_genome_bank and not os.path.isdir(template_genome_bank):
        logging.warning("Template genome bank directory not found: %s", template_genome_bank)
        sys.exit(1)

    if template_rerank_topn is not None and int(template_rerank_topn) < 0:
        logging.warning("Template recommendation requires '--template-rerank-topn' to be 0 or greater")
        sys.exit(1)

    validate_template_score_options(run_ns)

    if template_recommendation_only:
        if not getattr(run_ns, 'auto_template', False):
            logging.warning("Template recommendation-only mode requires '--auto-template'")
            sys.exit(1)
        if not pmr_generation:
            logging.warning("Template recommendation-only mode requires primary modeling option ('-p')")
            sys.exit(1)
        if smr_generation:
            logging.warning("Template recommendation-only mode cannot be combined with secondary modeling ('-s')")
            sys.exit(1)


def validate_template_score_options(run_ns):
    for field_name, default in _TEMPLATE_SCORE_DEFAULTS.items():
        value = float(getattr(run_ns, field_name, default))
        if value < 0.0 or value > 1.0:
            logging.warning(
                "Template recommendation weight '%s' must be within [0, 1]",
                _format_cli_option(field_name),
            )
            sys.exit(1)

    for left_field, right_field in _TEMPLATE_SCORE_PAIRS:
        left_value = float(getattr(run_ns, left_field, _TEMPLATE_SCORE_DEFAULTS[left_field]))
        right_value = float(getattr(run_ns, right_field, _TEMPLATE_SCORE_DEFAULTS[right_field]))
        if not math.isclose(left_value + right_value, 1.0, abs_tol=1e-9):
            logging.warning(
                "Template recommendation weights '%s' and '%s' must sum to 1.0",
                _format_cli_option(left_field),
                _format_cli_option(right_field),
            )
            sys.exit(1)


def _format_cli_option(field_name):
    return '--' + field_name.replace('_', '-')


def _read_executable_magic(candidate):
    try:
        with open(candidate, "rb") as handle:
            return handle.read(4)
    except OSError:
        return b""


def _is_unix_binary_compatible(candidate):
    magic = _read_executable_magic(candidate)
    if magic.startswith(b"#!"):
        return True
    if sys.platform == "darwin":
        return magic in _MACHO_MAGICS
    if sys.platform.startswith("linux"):
        return magic == _ELF_MAGIC
    return True


# Adopted from antismash.utils
def locate_executable(name):
    "Find an executable in the path and return the full path"
    valid_windows_suffixes = {".exe", ".bat", ".cmd"}

    def _is_executable(candidate):
        if not isfile(candidate):
            return False
        if sys.platform == 'win32':
            return os.path.splitext(candidate)[1].lower() in valid_windows_suffixes
        return os.access(candidate, os.X_OK) and _is_unix_binary_compatible(candidate)

    candidate_names = [name]
    if sys.platform == 'win32' and os.path.splitext(name)[1] == "":
        candidate_names = [name + ".exe", name + ".bat", name + ".cmd"]

    repo_root = abspath(join(dirname(__file__), os.pardir))
    search_paths = list(os.environ.get("PATH", "").split(os.pathsep))
    search_paths.extend([join(repo_root, "bin"), join(os.getcwd(), "bin")])

    for candidate_name in candidate_names:
        file_path, _ = split(candidate_name)
        if file_path != "":
            if _is_executable(candidate_name):
                logging.debug("Found executable %r", candidate_name)
                return candidate_name

        for p in search_paths:
            if not p:
                continue
            full_name = join(p, candidate_name)
            if _is_executable(full_name):
                logging.debug("Found executable %r", full_name)
                return full_name

    return None


# Adopted from antismash.utils
# Ignore the pylint warning about input being redifined, as we're just
# following the subprocess names here.
# pylint: disable=redefined-builtin
def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
         logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
         raise


# Adopted from antismash.utils
def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
         # force into a tuple
         types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features


# Adopted from antismash.utils
def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")


# Adopted from antismash.utils
def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"


# For regular update of the cache:
# KEGG updates its EC_number:reaction ID pairs time to time,
#and old caches can cause an error
def time_bomb(cache_file, config_ns):
    today = datetime.datetime.today()
    modified_date = datetime.datetime.fromtimestamp(getmtime(cache_file))
    file_age = today - modified_date

    if int(file_age.days) > int(config_ns.utils.time_bomb_duration):
        logging.debug('File %s is older than %s days (currently %s days)',
                cache_file, config_ns.utils.time_bomb_duration, file_age.days)
        os.remove(cache_file)
        logging.debug('File %s was removed', cache_file)
    else:
        logging.debug('File %s has not reached %s days (currently %s days)',
                cache_file, config_ns.utils.time_bomb_duration, file_age.days)


def get_keggid_from_mnxr(mnxr, io_ns, primary_model_ns):
    if len(io_ns.mnxr_kegg_dict[mnxr]) > 1:
        keggid_list = []

        for keggid in io_ns.mnxr_kegg_dict[mnxr]:
            if keggid in primary_model_ns.rxnid_info_dict:
                keggid_list.append(keggid)

        if len(keggid_list) == 1:
            kegg_id = keggid_list[0]
        # Choose KEGG reaction ID with a greater value for multiple KEGG IDs given to MNXR
        elif len(keggid_list) > 1:
            keggid_list.sort()
            kegg_id = keggid_list[-1]

    elif len(io_ns.mnxr_kegg_dict[mnxr]) == 1:
        kegg_id = io_ns.mnxr_kegg_dict[mnxr][0]

    return kegg_id


def check_duplicate_rxn(model, rxn2):

    duplicate_rxn = []

    for j in range(len(model.reactions)):
        rxn1 = model.reactions[j]

        comparison = compare_rxns(rxn1, rxn2)

        if comparison == 'same':
            duplicate_rxn.append(rxn1.id)

    if len(duplicate_rxn) >= 1:
        return 'duplicate'
    elif len(duplicate_rxn) == 0:
        return 'unique'


def compare_rxns(rxn1, rxn2):

    rxn1_metab_dict = {}
    for i in range(len(rxn1.metabolites)):
        rxn1_metab_dict[str(list(rxn1.metabolites.keys())[i])] = \
                    float(rxn1.metabolites[list(rxn1.metabolites.keys())[i]])

    rxn2_metab_dict = {}
    for i in range(len(rxn2.metabolites)):
        rxn2_metab_dict[str(list(rxn2.metabolites.keys())[i])] = \
                    float(rxn2.metabolites[list(rxn2.metabolites.keys())[i]])

    if rxn1_metab_dict.items() == rxn2_metab_dict.items():

        if rxn1.reversibility == rxn2.reversibility:
            return 'same'
        else:
            return 'different'
    else:
        return 'different'


#'add_reaction' requires writing/reloading of the model
def stabilize_model(model, folder, label, diff_name=False):

    if diff_name:
        model_name = '%s.xml' %label
    else:
        if label:
            model_name = 'model_%s.xml' %label
        else:
            model_name = 'model.xml'

    for rxn in model.reactions:
        if rxn.gene_reaction_rule == "()":
            rxn.gene_reaction_rule = ""

    ensure_modern_cobra_attrs(model)

    cobra.io.write_sbml_model(model, join('%s' %folder, model_name))
    model = cobra.io.read_sbml_model(join('%s' %folder, model_name))

    return model


#Output: a dictionary file for major Exchange reactions {Exchange reaction ID:flux value}
def get_exrxnid_flux(model, template_exrxnid_flux_dict):

    target_exrxnid_flux_dict = {}
    model.optimize()

    for exrxn_id in template_exrxnid_flux_dict:
        if exrxn_id in model.reactions:
            rxn = model.reactions.get_by_id(exrxn_id)
            target_exrxnid_flux_dict[exrxn_id] = rxn.flux
        else:
            continue
    return target_exrxnid_flux_dict


#Output: a list file having either T or F for major Exchange reactions
def check_exrxn_flux_direction(
        template_exrxnid_flux_dict, target_exrxnid_flux_dict, config_ns):

    exrxn_flux_change_list = []

    for exrxn_id in template_exrxnid_flux_dict:
        if exrxn_id in target_exrxnid_flux_dict:
            template_exrxn_flux = template_exrxnid_flux_dict[exrxn_id]
            target_exrxn_flux = target_exrxnid_flux_dict[exrxn_id]

            if float(template_exrxn_flux) == 0:
                if target_exrxn_flux == 0:
                    ratio_exrxn_flux = 1
                    logging.debug("%s: from zero to zero flux", exrxn_id)

                else:
                    ratio_exrxn_flux = False
                    logging.debug("%s: from zero to non-zero flux", exrxn_id)
                    continue

            else:
                ratio_exrxn_flux = float(target_exrxn_flux)/float(template_exrxn_flux)

            #Similar species are allowed to uptake nutrients within a decent range
            if ratio_exrxn_flux > 0 and ratio_exrxn_flux < float(config_ns.cobrapy.nutrient_uptake_rate) \
                and 1/ratio_exrxn_flux < float(config_ns.cobrapy.nutrient_uptake_rate):
                exrxn_flux_change_list.append('T')

            #Cause drastic changes in Exchange reaction fluxes
            #(direction and/or magnitude)
            else:
                logging.debug("Drastic change occured in %s: %f -> %f" \
                              %(exrxn_id, float(template_exrxn_flux), target_exrxn_flux))
                exrxn_flux_change_list.append('F')

    return exrxn_flux_change_list
