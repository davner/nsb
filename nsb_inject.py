#!/usr/bin/env python3
from astropy.io import ascii
from astropy.io import fits
# import tqdm for fancy progress bars!
from tqdm import tqdm
import numpy as np
import pandas as pd
import argparse
import pdb
def inject_headers_sqm(file, format):
    """
    Inject headers into the fits files with sqm measurements, telescope and
    zenith.

    Parameters
    file   | the sqm file (either nsb_sqm_tel.txt, nsb_sqm_zenith.txt)
    format | to know what header names to use
    """
    print('\n### Injecting sqm {} info into fits files'.format(format))
    data = ascii.read(file)

    if format == 'telescope':
        prefix = 'ST'
        do_az_elv = True
    elif format == 'zenith':
        prefix = 'SZ'
        do_az_elv = False

    # pass the data that will be iterated in the for loop in the tqdm function
    for row in tqdm(data):
        with fits.open(row['filename'], mode='update') as fits_file:
            fits_header = fits_file[0].header

            fits_header['{}-OBS'.format(prefix)] = (
                row['sqm_ut'], 'UTC SQM expsoure'
            )
            #fits_header['{}-TDEL'.format(prefix)] = (
                #row['t_delta'], 'time diff (seconds) fits UTC to SQM UTC'
            #)
            if do_az_elv:
                fits_header['{}-AZ'.format(prefix)] = (row['az'], 'degrees')
                fits_header['{}-ELV'.format(prefix)] = (row['elv'], 'degrees')

            fits_header['{}-COUNT'.format(prefix)] = (row['counts'], 'adu')
            fits_header['{}-NSB'.format(prefix)] = (
                row['reading'], 'nsb mag/arcsec^2'
            )

            fits_file.flush()

    print('Done...')

def inject_headers_boltwood(file):
    """
    Inject headers into the fits files with boltwood measurements.

    Parameters
    file | the boltwood file (nsb_bolt.txt)
    """
    print('\n### Injecting boltwood info into fits files')
    data = ascii.read(file)

    for row in tqdm(data):
        #pdb.set_trace()
        with fits.open(row['filename'], mode='update') as fits_file:
            fits_header = fits_file[0].header

            fits_header['BT-OBS'] = (
                row['bolt_ut'], 'UTC of boltwood measure'
            )
            #fits_header['BT-TDEL'] = (
                #row['t_del'], 'time diff (seconds) fits UTC to bolt UTC'
            #)
            fits_header['BT-SKYT'] = (row['sky_t'], 'sky temp (f)')
            fits_header['BT-AMBT'] = (row['amb_t'], 'ambient temp (f)')
            fits_header['BT-WIND'] = (row['wind'], 'wind (miles)')
            fits_header['BT-HUM'] = (row['hum'], 'humidity')
            fits_header['BT-DEW'] = (row['dew'], 'dew point')
            fits_header['BT-CLOUD'] = (
                row['cloud'], 'cloud flag 0=unk, 1=clr, 2=cloudy, 3=veryCloudy'
            )

            fits_file.flush()

    print('Done...')

def inject_headers_nsb(file):
    """
    Inject headers into the fits files with night sky brightness measurements.

    Parameters
    file | the nsb data file (nsb_data.txt)
    """
    print('\n### Injecting nsb measurements info fits files')
    data = ascii.read(file)

    for row in tqdm(data):
        with fits.open(row['filename'], mode='update') as fits_file:
            fits_header = fits_file[0].header

            fits_header['UA-MEAN'] = (
                row['mean_bkg'], 'mean background (counts/pixel)'
            )


            fits_header['UA-MED'] = (
                row['med_bkg'], 'median background (counts/pixel)'
            )
            fits_header['UA-STD'] = (row['std'], 'std of background')

            # check to see if value is nan, if so replace with 0
            if np.isnan(row['zp']) or np.isnan(row['nsb']):
                zp_value = (0, 'failed to find zeropoint')
                nsb_value = (0, 'failed to find nsb')
            else:
                zp_value = (row['zp'], 'zeropoint (mag)')
                nsb_value = (row['nsb'], 'nsb measurement (mag/arcsec^2)')

            fits_header['UA-ZP'] = zp_value
            fits_header['UA-NSB'] = nsb_value

            fits_file.flush()

    print('Done...')

def parse_arguments():
    """
    Create a parser to get the file from command line.

    Parameters
    ----------
    None

    Return
    ------
    args.files : the files in list form
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n',
        '--do_nsb',
        help='inject pomenis nsb data',
        action='store_true'
    )
    parser.add_argument(
        '-b',
        '--do_boltwood',
        action='store_true',
        help='inject boltwood data'
    )
    parser.add_argument(
        '-t',
        '--do_tel_sqm',
        action='store_true',
        help='inject telescope sqm readings'
    )
    parser.add_argument(
        '-z',
        '--do_zenith_sqm',
        action='store_true',
        help='inject zenith sqm readings'
    )

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    """
    Start of the program.
    """
    print('\n####### nsb_inject.py')

    args = parse_arguments()

    if args.do_tel_sqm:
        tel_sqm_file = '*sqm_tel*.csv'
        inject_headers_sqm(tel_sqm_file, 'telescope')

    if args.do_zenith_sqm:
        zenith_sqm_file = '*sqm_zenith*.csv'
        inject_headers_sqm(zenith_sqm_file, 'zenith')

    if args.do_boltwood:
        boltwood_file = '*bolt.csv'
        inject_headers_boltwood(boltwood_file)

    if args.do_nsb:
        nsb_file = '*nsb_data*.csv'
        inject_headers_nsb(nsb_file)

    print('\nFINISHED\n')
