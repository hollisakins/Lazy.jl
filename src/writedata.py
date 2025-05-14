from astropy.io import fits
import os

def write_data(filename, data, extname="DATA"):
    if not os.path.exists(filename):
        with fits.open(filename, mode='append') as hdul:
            # Create a new FITS file
            hdu = fits.PrimaryHDU()
            hdul.append(hdu)
            hdul.flush()
    
    with fits.open(filename, mode='append') as hdul:
        # Check if the extension already exists
        if extname in hdul:
            # If it exists, remove it
            hdul.remove(extname)
    
        if type(data) == dict:
            # Create a FITS table with the specified columns
            coldefs = []
            for col in data.keys():
                coldefs.append(
                    fits.Column(
                        name=col, 
                        format=data[col]["format"], 
                        unit=data[col]["unit"] if "unit" in data[col] else None, 
                        array=data[col]["data"]
                    )
                )
            # Create a new BinTableHDU with the specified data
            hdu = fits.BinTableHDU.from_columns(coldefs)
        else:
            # Create a new ImageHDU
            hdu = fits.ImageHDU(data=data, name=extname)

        hdul.append(hdu)
        hdul.flush()

    return