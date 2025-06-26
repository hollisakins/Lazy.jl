from astropy.io import fits
import os

def write_data(filename, data, extname):
    # If this is the first extension being written (SUMMARY), remove existing file
    if extname == "SUMMARY" and os.path.exists(filename):
        os.remove(filename)
    
    # Create file if it doesn't exist
    if not os.path.exists(filename):
        # Create a new FITS file with just a primary HDU
        hdu = fits.PrimaryHDU()
        hdu.writeto(filename)
    
    # Append the new extension
    with fits.open(filename, mode='append') as hdul:
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
            hdu = fits.BinTableHDU.from_columns(coldefs, name=extname)
        else:
            # Create a new ImageHDU
            hdu = fits.ImageHDU(data=data, name=extname)

        hdul.append(hdu)
        hdul.flush()

    return