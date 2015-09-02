

c_types = {
    'string' : 'string',
    'byte' : 'unsigned char',
    'integer*1' : 'char',
    'short' : 'short',
    'integer*2' : 'short',
    'integer' : 'int',
    'integer*4' : 'int',
    'long' : 'long',
    'integer*8' : 'long',
    'unsigned byte' : 'unsigned byte',
    'unsigned integer*1' : 'unsigned byte',
    'unsigned short' : 'unsigned short',
    'unsigned integer*2' : 'unsigned short',
    'unsigned integer' : 'unsigned integer',
    'unsigned integer*4' : 'unsigned integer',
    'unsigned long' : 'unsigned long',
    'unsigned integer*8' : 'unsigned long',
    'float' : 'float',
    'real' : 'float',
    'real*4' : 'float',
    'unsigned float' : 'unsigned float',
    'unsigned real' : 'unsigned float',
    'unsigned real*4' : 'unsigned float',
    'double' : 'double',
    'real*8' : 'double',
    'unsigned double' : 'unsigned double',
    'unsigned real*8' : 'unsigned double',
    'complex' : 'complex',
    'complex*8' : 'complex',
    'double complex' : 'double complex',
    'complex*16' : 'double complex'
}

fortran_types = {
    'string' : 'string',
    'byte' : 'integer*1',
    'integer*1' : 'integer*1',
    'short' : 'integer*2',
    'integer*2' : 'integer*2',
    'integer' : 'integer*4',
    'integer*4' : 'integer*4',
    'long' : 'integer*8',
    'integer*8' : 'integer*8',
    'unsigned byte' : 'unsigned integer*1',
    'unsigned integer*1' : 'unsigned integer*1',
    'unsigned short' : 'unsigned integer*2',
    'unsigned integer*2' : 'unsigned integer*2',
    'unsigned integer' : 'unsigned integer*4',
    'unsigned integer*4' : 'unsigned integer*4',
    'unsigned long' : 'unsigned integer*8',
    'unsigned integer*8' : 'unsigned integer*8',
    'float' : 'real*4',
    'real' : 'real*4',
    'real*4' : 'real*4',
    'unsigned float' : 'unsigned real*4',
    'unsigned real' : 'unsigned real*4',
    'unsigned real*4' : 'unsigned real*4',
    'double' : 'real*8',
    'real*8' : 'real*8',
    'unsigned double' : 'unsigned real*8',
    'unsigned real*8' : 'unsigned real*8',
    'complex' : 'complex*8',
    'complex*8' : 'complex*8',
    'double complex' : 'complex*16',
    'complex*16' : 'complex*16'
}

type_sizes = {
    'string' : 1,
    'byte' : 1,
    'integer*1' : 1,
    'short' : 2,
    'integer*2' : 2,
    'integer' : 4,
    'integer*4' : 4,
    'long' : 8,
    'integer*8' : 8,
    'unsigned byte' : 1,
    'unsigned integer*1' : 1,
    'unsigned short' : 2,
    'unsigned integer*2' : 2,
    'unsigned integer' : 4,
    'unsigned integer*4' : 4,
    'unsigned long' : 8,
    'unsigned integer*8' : 8,
    'float' : 4,
    'real' : 4,
    'real*4' : 4,
    'unsigned float' : 4,
    'unsigned real' : 4,
    'unsigned real*4' : 4,
    'double' : 8,
    'real*8' : 8,
    'unsigned double' : 8,
    'unsigned real*8' : 8,
    'complex' : 8,
    'complex*8' : 8,
    'double complex' : 16,
    'complex*16' : 16
}

def get_c_type (parsed_type):
    return c_types [parsed_type]

def get_fortran_type (parsed_type):
    return fortran_types [parsed_type]

def get_size (parsed_type):
    return type_sizes [parsed_type]


