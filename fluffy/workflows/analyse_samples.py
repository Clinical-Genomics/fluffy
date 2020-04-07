"""Workflow for analysing all samples in a batch"""

from pathlib import Path
from typing import Iterator

from fluffy.workflows.align import align_individual
from fluffy.workflows.amycne import estimate_ffy
from fluffy.workflows.picard import picard_qc_workflow
from fluffy.workflows.preface import preface_predict_workflow
from fluffy.workflows.wisecondor import wisecondor_xtest_workflow


def analyse_workflow(
    samples: Iterator[dict],
    out_dir: Path,
    configs: dict,
    skip_preface: bool = False,
    dry_run: bool = False,
):
    """Run the wisecondor chromosome x analysis"""
    jobids = []
    for sample in samples:
        sample_id = sample["sample_id"]
        sample_outdir = out_dir / sample_id
        # This will fail if dir already exists
        sample_outdir.mkdir(parents=True)
        # Should aligned files go here?

        align_jobid = align_individual(
            configs=configs, sample=sample, out_dir=out_dir, dry_run=dry_run
        )

        ffy_jobid = estimate_ffy(
            configs=configs,
            out_dir=out_dir,
            sample_id=sample_id,
            align_jobid=align_jobid,
        )
        jobids.append(ffy_jobid)

        picard_jobid = picard_qc_workflow(
            configs=configs,
            out_dir=out_dir,
            sample_id=sample_id,
            align_jobid=align_jobid,
        )
        jobids.append(picard_jobid)

        wcx_test_jobid = wisecondor_xtest_workflow(
            configs=configs,
            out_dir=out_dir,
            sample_id=sample_id,
            align_jobid=align_jobid,
        )

        jobids.append(wcx_test_jobid)

        if not skip_preface:
            preface_predict_jobid = preface_predict_workflow(
                configs=configs,
                out_dir=out_dir,
                sample_id=sample_id,
                dependency=wcx_test_jobid,
            )
            jobids.append(preface_predict_jobid)

    run_summarise = summarise(config, args)
    summarise_batch = Slurm(
        "summarise_batch-{}".format(args.project.strip("/").split("/")[-1]),
        {
            "account": config["slurm"]["account"],
            "partition": "core",
            "time": config["slurm"]["time"],
        },
        log_dir="{}/logs".format(args.out),
        scripts_dir="{}/scripts".format(args.out),
    )
    summarise_batch.run(run_summarise, depends_on=jobids)