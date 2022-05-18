from astropy.io import fits
import glob
import os
import shutil


def sort_data_targname_filt_propid(input_dir, output_dir, file_type,
                                   targname_mappings=None):

    """ Files in `input_dir` of type `file_type` are sorted into subdirectories
        in `output_dir` by target name, filter, and proposal ID.

        The outermost subdirectory is target name, the second is filter and
        the innermost directory containing the files is proposal ID.

        Parameters
        ----------
        input_dir : str
            Full path to where files to be sorted are currently located.
        output_dir : str
            Directory to sort files into in subdirectories.
         file_type : str
            Three-letter fits file extention (e.g flt, flc...). If 'any',
            all fits files in `input_dir` will be sorted.
         targname_mappings : None or dict
            If targets may go by different names in various files, provide
            a dictionary containing what their name should be mapped to, and
            the corresponding name variations. For example:

            targname_mappings = {'G191B2B' : ['G191B2B'],
                     'GD153' : ['GD153', 'GD-153'],
                     'GRW70' : ['GRW+70D5824', 'GRW+70D']}

            If None, the each file will be sorted into a subdirectory based
            on what`targname` is in each file header.

         """

    input_dir = os.path.join(input_dir, '')
    output_dir = os.path.join(output_dir, '')

    if file_type == 'any':
        file_type = '*'

    for f in glob.glob(input_dir + '*{}.fits'.format(file_type)):
        print(f)
        hdr = fits.open(f)[0].header
        targname = hdr['targname']
        if targname_mappings:  # get true name
            for key, val in targname_mappings.items():
                if targname in val:
                    targname = key
        proposid = str(hdr['proposid'])
        filt = hdr['filter']

        output_dirr = os.path.join(output_dir, targname, filt, proposid, '')

        if not os.path.isdir(output_dirr):
            print('Making directory {}.'.format(output_dirr))
            os.makedirs(output_dirr)

        print('Moving {} to {}'.format(f, output_dirr+os.path.basename(f)))
        shutil.move(f, output_dirr + os.path.basename(f))


def sort_data_targname_filt(input_dir, output_dir, file_type,
                            targname_mappings=None):

    """ Files in `input_dir` of type `file_type` are sorted into subdirectories
        in `output_dir` by target name, and filter.

        The outer directory is target name, and the inner directory containing
        the files is filter.

         Parameters
         ----------
        input_dir : str
            Full path to where files to be sorted are currently located.
        output_dir : str
            Directory to sort files into in subdirectories.
        file_type : str
            Three-letter fits file extention (e.g flt, flc...). If 'any',
            all fits files in `input_dir` will be sorted.
        targname_mappings : None or dict
            If targets may go by different names in various files, provide
            a dictionary containing what their name should be mapped to, and
            the corresponding name variations. For example:

            targname_mappings = {'G191B2B' : ['G191B2B'],
                     'GD153' : ['GD153', 'GD-153'],
                     'GRW70' : ['GRW+70D5824', 'GRW+70D']}

            If None, the each file will be sorted into a subdirectory based
            on what`targname` is in each file header.

         """

    input_dir = os.path.join(input_dir, '')
    output_dir = os.path.join(output_dir, '')

    if file_type == 'any':
        file_type = '*'

    for f in glob.glob(input_dir + '*{}.fits'.format(file_type)):
        print(f)
        hdr = fits.open(f)[0].header
        targname = hdr['targname']
        if targname_mappings:  # get true name
            for key, val in targname_mappings.items():
                if targname in val:
                    targname = key
        filt = hdr['filter']

        output_dirr = os.path.join(output_dir, targname, filt, '')

        if not os.path.isdir(output_dirr):
            print('Making directory {}.'.format(output_dirr))
            os.makedirs(output_dirr)

        print('Moving {} to {}'.format(f, output_dirr+os.path.basename(f)))
        shutil.move(f, output_dirr + os.path.basename(f))
