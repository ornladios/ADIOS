# Identity plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_identity_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_identity_read.c)

# Zlib plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_zlib_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_zlib_read.c)

# Bzip2 plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_bzip2_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_bzip2_read.c)

# Szip plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_szip_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_szip_read.c)

# ISOBAR plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_isobar_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_isobar_read.c)

# APLOD plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_aplod_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_aplod_read.c)

# ALACRITY plugin:
set(transforms_write_method_SOURCES ${transforms_write_method_SOURCES} transforms/adios_transform_alacrity_write.c)
set(transforms_read_method_SOURCES ${transforms_read_method_SOURCES} transforms/adios_transform_alacrity_read.c)
