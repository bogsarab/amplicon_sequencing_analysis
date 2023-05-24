from pydantic.main import BaseModel
from pydantic import BaseSettings
from typing import List, Dict


class BarcodesConfig(BaseModel): # barcodes=run descr barcodes
    CGATGT: str
    CAGTAC: str
    TGACCA: str
    TTAGGC: str


class PrimersConfig(BaseModel):
    forward: str
    reverse: str


class FileConfig(BaseSettings):
    file_path: str
    working_dir: str
    out_dir: str
    figures_dir: str

    class Config:
        env_file = '.env'
        case_sensitive = False


class RunSteps(BaseSettings):
    run_preprocess: bool
    run_analysis: bool
    run_indel_analysis: bool

    class Config:
        env_file = '.env'
        case_sensitive = False


class RunConfig(BaseModel):

    run_steps: RunSteps
    barcodes: BarcodesConfig
    primers: PrimersConfig
    samples: List[str]
    ref_seq: str
    original_seq: str
    possible_bases: str
