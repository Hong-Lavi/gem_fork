import sys

import cobra

import run_gmsm


def test_primary_modeling_without_ec_uses_pruned_model(monkeypatch, tmp_test_dir):
    captured = {}
    pruned_model = cobra.Model("pruned_model")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_gmsm.py",
            "-i",
            "dummy_input.gbk",
            "-p",
            "-d",
            "-o",
            tmp_test_dir,
        ],
    )

    monkeypatch.setattr(run_gmsm, "make_folder", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm.utils, "setup_logging", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        run_gmsm,
        "setup_outputfolders",
        lambda run_ns, io_ns: setattr(io_ns, "outputfolder2", tmp_test_dir),
    )
    monkeypatch.setattr(run_gmsm.utils, "check_input_options", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm, "show_input_options", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm, "check_input_filetype", lambda *args, **kwargs: "genbank")
    monkeypatch.setattr(run_gmsm, "load_config", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm, "check_prereqs", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm, "remove_tmp_model_files", lambda *args, **kwargs: None)

    def fake_get_target_genome_from_input(filetype, run_ns, io_ns):
        io_ns.targetGenome_locusTag_aaSeq_dict = {"gene1": "M"}
        io_ns.targetGenome_locusTag_ec_dict = {}

    monkeypatch.setattr(run_gmsm, "get_target_genome_from_input", fake_get_target_genome_from_input)
    monkeypatch.setattr(run_gmsm, "get_fasta_files", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm, "get_homologs", lambda *args, **kwargs: None)
    monkeypatch.setattr(run_gmsm, "get_pickles_prunPhase", lambda *args, **kwargs: cobra.Model("template_model"))
    monkeypatch.setattr(run_gmsm, "run_prunPhase", lambda *args, **kwargs: pruned_model)
    monkeypatch.setattr(run_gmsm, "prune_unused_metabolites", lambda model: None)

    def fake_generate_outputs(folder, runtime, run_ns, io_ns, homology_ns, primary_model_ns, secondary_model_ns, **kwargs):
        captured["model"] = kwargs["cobra_model"]

    monkeypatch.setattr(run_gmsm, "generate_outputs", fake_generate_outputs)

    run_gmsm.main()

    assert captured["model"] is pruned_model
