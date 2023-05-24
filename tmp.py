from configs import FileConfig, RunConfig, BarcodesConfig
from pprint import pprint

fileconfig = FileConfig()
runconfig = RunConfig.parse_file('run_descriptor.json')


pprint(fileconfig.__dict__)
pprint(runconfig.samples)

