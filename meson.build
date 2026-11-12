project('kagome_AC_CMC', ['cpp', 'c'],
  version: '0.0.1',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=release', 'cpp_std=c++2a' ]
)


json_dep = dependency('nlohmann_json', required: true)

latlib_dep = dependency('lattice_lib', version: '>= 2.3.0')

fftw_dep = dependency('fftw3', required: true)

hdf5_dep = dependency('hdf5', required: true, language: 'cpp')

cc=meson.get_compiler('c')

# Detect whether HDF5 was built with MPI support
if cc.compiles('''
    #include <hdf5.h>
    int main() { return 0; }
''',
    name: 'HDF5 compiles without MPI',
    dependencies: hdf5_dep
)

    message('Serial HDF5 detected')
else
  # Parallel HDF5 pulls in MPI, so require it and fold it into the dependency.
    mpi_dep = dependency('mpi', include_type: 'system', language: 'cpp', required: true)
    hdf5_dep = declare_dependency(
        dependencies: [hdf5_dep, mpi_dep]
    )
    message('parallel HDF5 detected')
endif

CMC_inc=[
  include_directories(['include'])
  ]

g_deps = [
    latlib_dep,
    json_dep,
    fftw_dep,
    hdf5_dep
    ]

mc_lib = static_library('cmc_lib',
  files('src/MC.cpp'),
  dependencies: g_deps,
  include_directories: CMC_inc
  )

anneal_bin = executable('anneal', 
  files('src/anneal.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: CMC_inc,
  install: true
  )

AC_sus_bin = executable('AC_sus', files('src/AC_sus.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: CMC_inc,
  install: true
  )

DC_bin = executable('DC', files('src/DC.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: CMC_inc,
  install: true
  )

eigen_dep = dependency('eigen3', required: true, include_type: 'system')

ltgs_bin = executable('ltgs',
  files('src/luttinger_tisza.cpp'),
  dependencies: g_deps + [eigen_dep],
  link_with: mc_lib,
  include_directories: CMC_inc,
  install: true
  )

py_scripts = ['ltgs_plot', 'ltgs_sweep']
script_dir = get_option('datadir') / 'kagome_cmc' / 'scripts'

if get_option('dev')
  # Point directly at the source tree so edits take effect without reinstalling.
  foreach s : py_scripts
    install_symlink(s,
      pointing_to: meson.project_source_root() / 'results' / s + '.py',
      install_dir: get_option('bindir'),
    )
  endforeach
else
  install_data(
    files(
      'results/ltgs_plot.py',
      'results/ltgs_sweep.py',
    ),
    install_dir: script_dir,
    install_mode: 'rwxr-xr-x',
  )
  # Absolute, prefix-derived target so the symlink resolves under whatever
  # prefix/datadir was configured (join drops prefix if datadir is absolute).
  foreach s : py_scripts
    install_symlink(s,
      pointing_to: get_option('prefix') / script_dir / s + '.py',
      install_dir: get_option('bindir'),
    )
  endforeach
endif


# subdir('test')

