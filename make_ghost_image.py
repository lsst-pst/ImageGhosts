# Generate an image that includes the ghosts from bright stars.

import os
import math
import numpy
import argparse
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.math as afwMath
from lsst.skymap.detail import WcsFactory

### Notes:
# If you're not actually using the WCS to get x/y locations for the stars (which didn't work for me),
# then the RA/Dec of the center of the WCS is irrelevant as long as it's consistent. Similar for the pixScale.
#
# And while I've just defined the RaCen/DecCen of the pointing (which is used to determine the distance of the star
# from the center of the fov and the rotation angle) as the mean of the input stars, you could do this differently.


### Dealing with ghost and image data.

def readGhost(ghostfilename, ghostDir, removeSource=True):
    # Read the text file format ghost files.
    usecols = numpy.arange(1, 1297, 1, 'int')
    # Read numpy array data.
    ghost = numpy.loadtxt(os.path.join(ghostDir, ghostfilename), skiprows=24, usecols=usecols, dtype=numpy.float32)
    # Note that the raw ghost images have a ~4pix 'point' that corresponds to the star.
    #  This is direct light and not part of the ghost.
    # Remove this direct light.
    if removeSource:
        # This seems like a good rule of thumb .. the direct light >.1, ghost <<1 (from looking at files). (and I tested it - this removes 4 pixels for our sample)
        ghost = numpy.where(ghost>0.2, 0, ghost)
    # Turn into an LSST afwImage.
    ghostim = afwImage.ImageF(ghost)
    return ghostim

def readAllGhosts(ghostdir='ghost_data', removeSource=True):
    # Look for the ghost files in the 'ghostdir' directory.
    tmp = os.listdir(ghostdir)
    ghostfiles = []
    d_angs = []
    for t in tmp:
        if ((t.startswith('LSST_Ghost')) & (t.endswith('.txt'))):
            ghostfiles.append(t)
            d_angs.append(t.split('_')[2].strip('.txt'))
    d_angs = numpy.array(d_angs, 'float')
    # d_angs in the filenames are just the steps across the focal plane.
    # Translate from steps into angles measured from the center.
    d_angs = d_angs * 0.02
    ghostIms = {}
    for ghostfile, d_ang in zip(ghostfiles, d_angs):
        ghostIms[d_ang] = readGhost(ghostfile, ghostdir, removeSource=removeSource)
    return ghostIms, d_angs

def makeWcs(pixScale, nx, ny, rotation=0., raCen=180.0, decCen=0.0):
    """Create a TAN WCS with pixScale/nx/ny, plus rotation (degrees) and RA/Dec center (degrees). """
    #Hardcoding TAN projection since this application is for focal plane size at most
    wcsFactory = WcsFactory(pixScale, 'TAN', rotation=rotation*afwGeom.radians)
    crpix_x = float(nx-1)/2.
    crpix_y = float(ny-1)/2.
    crval = afwCoord.Coord(afwGeom.Angle(raCen*afwGeom.degrees), afwGeom.Angle(decCen*afwGeom.degrees))
    return wcsFactory.makeWcs(afwGeom.Point2D(crpix_x, crpix_y), crval)

def makeDest(nx, ny, destPixScale, raCen, decCen):
    #Grow the output image so the input image fits no matter what the rotation
    #destSize = int(2*math.hypot(nx/2., ny/2.))
    #destBbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(destSize, destSize))
    # Keep the original size image.
    destBbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(nx, ny))
    destImage = afwImage.ImageF(destBbox)
    #Make WCS for final (possibly grown) image
    destWcs = makeWcs(destPixScale, destBbox.getDimensions().getX(), destBbox.getDimensions().getY(),
                      raCen=raCen, decCen=decCen)
    return destImage, destWcs

### Deal with star data and locations.

def readStars(starfile):
    # Read data.
    # Star's RA and Dec values are in degrees.
    file = open(starfile, 'r')
    stars = {}
    stars['ra'] = []
    stars['dec'] = []
    stars['mag'] = []
    for line in file:
        try:
            rai, deci, magi = line.split()
        except:
            print 'could not read', line
            continue
        stars['ra'].append(rai)
        stars['dec'].append(deci)
        stars['mag'].append(magi)
    file.close()
    for k in stars.keys():
        stars[k] = numpy.array(stars[k], 'float')
    # Calculate projection. Assume boresight is ~center of file.
    stars['raCen'] = stars['ra'].mean()
    stars['decCen'] = stars['dec'].mean()
    return stars

def calcStarPos(stars):
    deg2rad = numpy.pi/180.0
    ## Use spherical trig to just calculate d_ang and rotAng directly.
    dlon = (stars['dec'] - stars['decCen'])*deg2rad
    stars['d_ang'] = numpy.arccos(numpy.sin(stars['ra']*deg2rad)*numpy.sin(stars['raCen']*deg2rad) + \
                                  (numpy.cos(stars['ra']*deg2rad)*numpy.cos(stars['raCen']*deg2rad)*numpy.cos(dlon)))
    stars['d_ang'] /= deg2rad
    stars['rotAng'] = numpy.arctan2((numpy.sin(dlon)*numpy.cos(stars['ra']*deg2rad)),
                                    (numpy.cos(stars['ra']*deg2rad)*numpy.sin(stars['raCen']*deg2rad) -
                                     numpy.sin(stars['ra']*deg2rad)*numpy.cos(stars['raCen']*deg2rad)*numpy.cos(dlon)))
    stars['rotAng'] /= deg2rad
    return stars

if __name__ == '__main__':
    # Parse arguments from command line.
    parser = argparse.ArgumentParser(description="Create an 'image' consisting of ghost images resulting from stars.")
    parser.add_argument('--outImage', action="store", type=str, default='testOut.fits',
                   help='Name of the output (fits) image. (default=testOut.fits).')
    parser.add_argument('--ghostDir', action="store", type=str, default='.',
                     help='Location of the (input) ghost image (txt files) directory (default=.)')
    parser.add_argument('--starCat', action="store", type=str, default='star_cats/Monet_102_0.dat',
                     help='Name of the star catalog input file (default=star_cats/Monet_102_0.dat).')
    parser.add_argument('--nStars', action="store", type=int, default=-1,
                        help='Number of stars (sorted by magnitude) to use from catalog (<0 for all) (default=-1).')
    args = parser.parse_args()

    print 'Starting to read ghost files.'
    # Read ghost images & convert to afwImages.
    ghostdir = '.'
    ghostIms, d_angs = readAllGhosts(args.ghostDir)
    #print d_angs
    d_angs = numpy.sort(d_angs)
    # Calculate largest angle from center that should use these ghosts.
    d_ang_max_use = d_angs[-1:] + (d_angs[-1:] - d_angs[-2:-1])[0]
    print 'Read all ghost files.'

    # Write out ghosts in fits file for examination.
    #for d_ang in d_angs:
    #    ghostIm = ghostIms[d_ang]
    #    ghostIm.writeFits('Ghost_%.2f.fits' %(d_ang))

    # Read star file.
    stars = readStars(args.starCat)

    ## Set up final image.
    # Use one of the input images as a base for dimensions.
    d_ang = d_angs[0]
    nx, ny = ghostIms[d_ang].getDimensions()
    # pixel scale in ghost images is 10" per pixel
    pixScale = afwGeom.Angle(10.0/60./60./1.)
    # Create the destination image.
    destIm, destWcs = makeDest(nx, ny, pixScale, stars['raCen'], stars['decCen'])
    destArr = destIm.getArray()

    # Calculate star positions.
    stars = calcStarPos(stars)
    #print stars['d_ang'].min(), stars['d_ang'].max(), stars['rotAng'].min(), stars['rotAng'].max()

    # Set up warper to map input image to output WCS.
    # There are other warping kernels, but lanczos2 is the standard one
    warper = afwMath.Warper("lanczos2")

    idx = numpy.argsort(stars['mag'])
    if args.nStars > 0:
        idx = idx[0:args.nStars]
    for i in (idx):
        ## For each star, choose appropriate ghost & rotate. Add.
        # Find the closest ghost image.
        d_ang = stars['d_ang'][i]
        # Skip star if too far outside the last ghost image.
        if d_ang > d_ang_max_use:
            continue
        d = abs(d_angs - d_ang)
        didx = numpy.where(d == d.min())
        d_ang_closest = d_angs[didx][0]
        # Calculate scale for output image (assuming ghost source = 0 mag)
        fscale = numpy.power(10.0, -0.4*stars['mag'][i])
        # Change input WCS (particularly rotation angle), to match location of star.
        srcWcs = makeWcs(pixScale, nx, ny, rotation=stars['rotAng'][i],
                         raCen=stars['raCen'], decCen=stars['decCen'])
        # Warp ghost image to rotation angle.
        outIm = warper.warpImage(destWcs, ghostIms[d_ang_closest], srcWcs,
                                 destBBox=destIm.getBBox(afwImage.PARENT))
        #Coadding on array since we need to map the nan pixels to zero because
        #nan + 1 = nan
        destArr += numpy.nan_to_num(outIm.getArray()*fscale)
        print 'Added ghost for star %d (mag %f) with d_ang %.2f and rotAng %f' %(i, stars['mag'][i],
                                                                               d_ang_closest, stars['rotAng'][i])

    print 'Destination dimensions', destIm.getDimensions()
    destIm.writeFits(args.outImage)

