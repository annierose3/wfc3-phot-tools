"""
This module contains functions to search for and
download WFC3 data from MAST using Astroquery.

Authors
-------
    Mariarosa Marinelli, 2022-2023
    Clare Shanahan, Oct 2019

Notes
-----
    Data should be public access for calibration, but
    proprietary data will need access authentication.

Functions
---------
download_products(query_products, output_dir)
query_by_data_id(dataset_ids, file_type)
query_by_propid_targ_filter(prop_ids, target_names, filters, file_types)

"""
import os
import glob
import shutil
from astroquery.mast import Observations
from astropy.table import Table, vstack


def download_products(query_products, output_dir=os.getcwd()):
    """Download data products to `output_dir`.

    Parameters
    ----------
    query_products : `astropy.table.Table`
        Table of data products to download.
    output_dir : str
        Output directory where data products should be
        downloaded. If `output_dir` does not exist or is
        not specified, this function will default to the
        current working directory.

    Notes
    -----
    Files are initially downloaded to a temporary directory
    within `output_dir` called 'temp', so if a subdirectory
    'temp' already exists within `output_dir`, an error is
    raised.

    Usage
    -----
    Implemented in WFC3/UVIS staring mode pipeline.
    """
    if not os.path.exists(output_dir):
        output_dir = os.getcwd()

    temp_dir = os.path.join(output_dir, 'temp')

    if os.path.exists(temp_dir):
        raise Exception("'temp' directory already exists in `output_dir`: "\
                        f"{output_dir}")

    else:
        os.makedirs(temp_dir)
        print(f'Created temporary directory at: {temp_dir}')

    print(f'Downloading {len(query_products)} files.')
    Observations.download_products(query_products,
                                   download_dir=temp_dir,
                                   mrp_only=False)

    # move files from MAST download directories in `temp` to `output_dir`
    files = glob.glob(output_dir + 'temp/mastDownload/HST/*/*.fits')
    if len(files) > 0:
        print("Cleaning up 'temp' directory.")
        for f in files:
            shutil.copy(f, output_dir+os.path.basename(f))

        # remove temp directory
        shutil.rmtree(temp_dir)


def query_by_propid_targ_filter(prop_ids, target_names='any', filters='any',
                                file_types='any'):
    """Query MAST by program, target, and filter.

    This function uses Astroquery to query MAST for HST
    data from one or more specified program. Additional
    options include filtering on target names, detector
    filters, and file types. Returns a table of details for
    data products that are available to download through
    MAST.

    Parameters
    ----------
    prop_ids : str or list of str
        Proposal ID(s)
    target_names : str list of str
        Exact target names that should be returned in
        query, others will be excluded. All spelling/name
        variations that might appear in MAST should be
        provided. If 'any', all available targets will be
        returned.
    filters: str or list of str
        Filters that should be returned in query; others
        will be excluded. If 'any', all available filters
        will be returned.
    file_types : str or list of str
        File extention type(s) desired (i.e flt, flc,
        drz...), as a string for a single type or a list
        for many. If 'any', all available file types will
        be returned.

    Returns
    --------
    query_products : `astropy.table.Table`
        Table of products returned from query.

    Usage
    -----
    Implemented in WFC3/UVIS staring mode pipeline.
    """
    if type(prop_ids) != list:
        prop_ids = [prop_ids]
    if target_names == 'any':
        target_names = '*'
    if filters == 'any':
        filters = '*'

    query_products_total = Table()
    j = 0
    for i, prop_id in enumerate(prop_ids):  # iterate to avoid timeout
        print('Querying for data from {}.'.format(prop_id))

        obsTable = Observations.query_criteria(obs_collection='HST',
                                               proposal_id=prop_id,
                                               target_name=target_names,
                                               filters=filters)
        if file_types == 'any':
            query_products = Observations.get_product_list(obsTable)
        else:
            if type(file_types) == str:
                file_types = file_types.upper()
            if type(file_types) == list:
                file_types = [x.upper() for x in file_types]
            query_products = Observations.get_product_list(obsTable)
            query_products = Observations.filter_products(query_products,
                                                          productSubGroupDescription=file_types)

        if len(query_products) == 0:
            print('No records found in query.')
            j = 0
            if len(prop_ids) == 1:
                return
            continue

        if (i == 0) & (j == 0):
            query_products_total = query_products
        else:
            query_products_total = vstack([query_products,
                                           query_products_total])
        print('{} records found'.format(len(query_products)))
        j = 1

    return query_products_total


def query_by_data_id(dataset_ids, file_type):
    """Query MAST by file rootnames/ASN ID.

    This function uses Astroquery to query MAST for HST
    data matching one or more rootnames/ASN IDs. As an
    additional option, file types can be specified to limit
    what is returned by the query. Returns a table of
    details for data products available to download through
    MAST that match IDs in `dataset_ids`.

    Parameters
    ----------
    dataset_ids : str or list of str
        The dataset ID or IDs for which to query.
    file_type : str or list of str
        Type of files (ex. 'FLT', 'DRC') to be returned by
        the MAST query.

    Returns
    -------
    query_products_total : `Astropy.table.Table`
        Table of all data products matching given query
        parameters.

    Notes
    -----
    Does not appear to be implemented in any of the WFC3
    photometry monitoring pipelines. Deprecate?
    """
    # Initial query for all files in visit, since you can
    # only query by ASN ID and `dataset_ids` might contain
    # single exposure rootnames
    if type(dataset_ids) == str:
        visit_ids = [dataset_ids[0:6] + '*']
    if type(dataset_ids) == list:
        visit_ids = list(set([x[0:6]+'*' for x in dataset_ids]))

    query_products_total = None
    for i, idd in enumerate(visit_ids):
        obsTable = Observations.query_criteria(obstype='all',
                                               obs_collection='HST',
                                               obs_id=idd)
        if file_type == 'any':
            query_products = Observations.get_product_list(obsTable)
        else:
            if type(file_type) == str:
                file_type = file_type.upper()
            if type(file_type) == list:
                file_type = [x.upper() for x in file_type]
            query_products = Observations.get_product_list(obsTable)
            query_products = Observations.filter_products(query_products,
                                                          productSubGroupDescription=file_type)

        if i == 0:
            query_products_total = query_products
        else:
            query_products_total = vstack([query_products,
                                           query_products_total])

    # Initially all visit files were returned. Now select
    # only specified IDs.
    remove_rows = []
    for i, obs_id in enumerate(query_products_total['obs_id']):
        if obs_id not in dataset_ids:
            remove_rows.append(i)
    query_products_total.remove_rows(remove_rows)

    print('{} records found'.format(len(query_products_total)))
    return query_products_total
