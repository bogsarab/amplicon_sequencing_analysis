from configs import RunConfig, FileConfig
from preprocessing import preprocessing
from analysis import analysis
from indels import indel_analysis
import os


def main():
    os.chdir('C:/Users/Bognár Sarolta/Documents/labor/processzált/amp seq analysis/amplicon_seq_analysis-main')
    runconfig = RunConfig.parse_file("run_descriptor.json")
    fileconfig = FileConfig()

    if runconfig.run_steps.run_preprocess is True:
        preprocessing(runconfig=runconfig,   # pl itt állítottuk be h csak a preprocesset futtassa le
                      fileconfig=fileconfig)

    if runconfig.run_steps.run_analysis is True:
        analysis(runconfig=runconfig,
                 fileconfig=fileconfig)

    if runconfig.run_steps.run_indel_analysis is True:  
        indel_analysis(runconfig=runconfig,
                       fileconfig=fileconfig)


if __name__ == '__main__': 
    main()
