project(
  'stringspinner',
  'cpp', 'fortran',
  meson_version: '>=1.2',
  default_options: {
    'cpp_std':   'c++11',
    'buildtype': 'release',
    'libdir':    'lib',
  },
  version: 'Comput.Phys.Commun. 292 (2023) 108886', # no version number, so use the publication reference instead; arXiv: https://arxiv.org/abs/2305.05058
)

# dependency pythia8
pythia_link_args = []
pythia_config_prog = find_program('pythia8-config', required: false)
if not pythia_config_prog.found()
  error('Pythia8 is not found')
endif
foreach arg : run_command(pythia_config_prog, '--libs', check: true).stdout().strip().split()
  if not arg.contains('-rpath')
    pythia_link_args += arg
  endif
endforeach
pythia_version_res = run_command(pythia_config_prog, '--version', check: false)
pythia_version = pythia_version_res.stdout().strip()
if pythia_version_res.returncode() != 0 or pythia_version == ''
  warning('pythia version auto-detection failed, since it is likely older than 8.312; setting version to 8.3')
  pythia_version = '8.3' # fallback for versions < 8.312; see https://gitlab.com/Pythia8/releases/-/issues/448
endif
stringspinner_pythia_dep = declare_dependency(
  include_directories: run_command(pythia_config_prog, '--includedir', check: true).stdout().strip(),
  link_args: pythia_link_args,
  version: pythia_version,
)

# headers
