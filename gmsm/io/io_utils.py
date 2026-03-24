
import logging
import os
import sys


def _resolve_feature_gene_id(feature):
    for qualifier_name in ('locus_tag', 'protein_id', 'gene'):
        values = feature.qualifiers.get(qualifier_name)
        if values:
            return values[0].replace("-", "_")
    return None


def _is_bgc_feature(feature):
    return feature.type in ('region', 'cluster')


def _normalize_path(path):
    return path.replace(os.sep, '/')


def should_ignore_input_gbk_ec_annotations(run_ns):
    return bool(
        getattr(run_ns, 'ec_file', False)
    )


def get_antismash_version_from_gbk(seq_record, io_ns):
    anti_data = seq_record.annotations.get('structured_comment', {}).get('antiSMASH-Data', {})
    version = anti_data.get('Version')
    if version:
        io_ns.anti_version = int(str(version).split('.')[0])
        return

    bgc_feature_types = {feature.type for feature in seq_record.features}
    if 'cluster' in bgc_feature_types:
        io_ns.anti_version = 4
    elif 'region' in bgc_feature_types:
        io_ns.anti_version = 5
    else:
        io_ns.anti_version = None


def get_features_from_gbk(seq_record, run_ns, io_ns):

    seq_record_BGC_num_list = []
    seq_record_BGC_num_list.append(seq_record)
    BGC_num = 0

    ignore_embedded_ec = should_ignore_input_gbk_ec_annotations(run_ns)

    for feature in seq_record.features:
        if feature.type == 'CDS':

            # Prefer locus_tag, but accept protein_id or gene for GenBank inputs
            # that omit locus_tag qualifiers entirely.
            locusTag = _resolve_feature_gene_id(feature)
            if locusTag is None:
                logging.error("No usable gene identifier found in gbk file")
                sys.exit(1)

            # Note that the numbers of CDS and "translation" do not match.
            # Some CDSs do not have "translation".
            if feature.qualifiers.get('translation'):
                translation = feature.qualifiers.get('translation')[0]
                io_ns.targetGenome_locusTag_aaSeq_dict[locusTag] = translation

            # Used to find "and" relationship in the GPR association
            if feature.qualifiers.get('product'):
                # It is confirmed that each locus_tag has a single '/product' annotation.
                # Thus, it's OK to use '[0]'.
                product = feature.qualifiers.get('product')[0]
                io_ns.targetGenome_locusTag_prod_dict[locusTag] = product
                
            if not ignore_embedded_ec:
                if feature.qualifiers.get('EC_number'):
                    # Multiple 'EC_number's may exit for a single CDS.
                    # Never use '[0]' for the 'qualifiers.get' list.
                    ecnum = feature.qualifiers.get('EC_number')
                    io_ns.targetGenome_locusTag_ec_dict[locusTag] = ecnum

        if _is_bgc_feature(feature):
            io_ns.total_region += 1
            BGC_num += 1

    seq_record_BGC_num_list.append(BGC_num)
    io_ns.seq_record_BGC_num_lists.append(seq_record_BGC_num_list)


def get_features_from_fasta(seq_record, io_ns):
    locusTag = seq_record.id
    io_ns.targetGenome_locusTag_aaSeq_dict[locusTag] = seq_record.seq
    io_ns.targetGenome_locusTag_prod_dict[locusTag] = seq_record.description


def get_target_fasta(io_ns):

    if io_ns.targetGenome_locusTag_aaSeq_dict:
        target_fasta_dir = os.path.join(
                io_ns.outputfolder1, 'targetGenome_locusTag_aaSeq.fa')
        with open(target_fasta_dir,'w') as f:
            for locusTag in io_ns.targetGenome_locusTag_aaSeq_dict.keys():
                print('>%s\n%s' \
                %(str(locusTag), str(io_ns.targetGenome_locusTag_aaSeq_dict[locusTag])), file=f)
        io_ns.target_fasta = _normalize_path(target_fasta_dir)
    else:
        logging.warning("FASTA file 'targetGenome_locusTag_aaSeq.fa' not found")


#Look for pre-stored fasta file of the template model
def get_temp_fasta(run_ns, io_ns):

    for root, _, files in os.walk('./gmsm/io/data/input1/%s/' %(run_ns.orgName)):
        for f in files:
            if f.endswith('.fa'):
                tempFasta = os.path.join(root, f)
                io_ns.input1 = _normalize_path(root).rstrip('/') + '/'
                io_ns.temp_fasta = _normalize_path(tempFasta)

    if io_ns.temp_fasta:
        logging.debug("FASTA file for '%s' found", run_ns.orgName)
    else:
        logging.warning("FASTA file for '%s' not found", run_ns.orgName)
