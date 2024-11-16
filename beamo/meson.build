project('pythia_test', 'cpp')

pythia_config = find_program('pythia8-config')
pythia_link_args = []
foreach arg : run_command(pythia_config, '--libs', check: true).stdout().strip().split()
  if not arg.contains('-rpath')
    pythia_link_args += arg
  endif
endforeach
pythia_dep = declare_dependency(
  include_directories: run_command(pythia_config, '--includedir', check: true).stdout().strip(),
  link_args: pythia_link_args,
)

exe = executable('pythia_test', 'test.cpp', dependencies: [pythia_dep])
test('pythia_test', exe)
