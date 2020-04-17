"""An API to slurm"""

import logging
from pathlib import Path

from slurmpy import Slurm

LOG = logging.getLogger(__name__)


class SlurmAPI:
    """An API to SLURM"""

    def __init__(self,slurm_settings: dict, out_dir: Path):
        super(SlurmAPI, self).__init__()
        LOG.info("Initializing a slurm API")
        account=None
        if "account" in slurm_settings:
            account=slurm_settings["account"]
        if not isinstance(account, str):
            raise SyntaxError("Invalid account {}".format(account))
        self.account = account

        time=None
        if "time" in slurm_settings:
            time=slurm_settings["time"]

        if not isinstance(time, str):
            raise SyntaxError("Invalid time {}".format(time))
        if not len(time.split(":")) == 3:
            LOG.warning("Specify time on format hh:mm:ss")
            raise SyntaxError("Invalid time format {}".format(time))
        self.time = time
        if not isinstance(out_dir, Path):
            raise SyntaxError("Invalid out dir format {}".format(out_dir))

        self.log_dir = out_dir / "logs"
        self.scripts_dir = out_dir / "scripts"
        self.slurm_settings=slurm_settings
        self.job = None

    def create_job(self, name: str) -> Slurm:
        """Create a job for submitting to SLURM"""
        LOG.info("Create a slurm job with name %s", name)
        job = Slurm(
            name,
            self.slurm_settings,
            log_dir=str(self.log_dir),
            scripts_dir=str(self.scripts_dir),
        )
        return job

    def run_job(
        self, name: str, command: str, dependencies: list = None, dry_run: bool = False
    ) -> int:
        """Create and submit a job to slurm"""
        job = self.create_job(name=name)
        LOG.info("Submitting commands %s", command)
        if dependencies:
            LOG.info("Adding dependencies: %s", ",".join( str(dependency) for dependency in dependencies))
        jobid = 1
        if not dry_run:
            jobid = job.run(command, depends_on=dependencies)
        LOG.info("Submitted job %s with job id: %s", name, jobid)
        return jobid

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(account={self.account!r}, time={self.time!r}, "
            f"log_dir={self.log_dir!r}, scripts_dir={self.scripts_dir!r})"
        )
