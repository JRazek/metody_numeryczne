project(
    'jr_numeric', 'c', 'cpp', 
    version: '0.1.0',
    default_options : ['c_std=c18', 'cpp_std=c++20']
)

cpp = meson.get_compiler('cpp')

cpp_build_args=[
]
cpp_link_args=[
]

incdir = include_directories(
  'include/',
)

if get_option('build')=='release'
    cpp_build_args += ['-O3', '-Wall', '-DNDEBUG', '-g0']
else
	cpp_build_args += ['-g3', '-Wall', '-O0', '-fsanitize=address']
	cpp_link_args += ['-fsanitize=address', '-flto']
endif

fmt_dep = dependency('fmt', version: '>=9.0.0')

numeric_lib_dep = declare_dependency(include_directories: incdir, dependencies: [fmt_dep])

if host_machine.system() == 'linux'
  install_script = join_paths(meson.current_source_dir(), 'meson/install.sh')
  include_path = join_paths(meson.current_source_dir(), 'include')

  install_path = '/usr/local/include/jr_numeric'

  meson.add_install_script(install_script, include_path, install_path)
endif


if get_option('build_examples')==true
  subdir('example')
endif
