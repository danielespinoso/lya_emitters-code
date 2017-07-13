import numpy as np
import webbrowser as wb

def browse_jplus_images(data):
  BrowseObjImages = True
  if BrowseObjImages:
    fracbrowse = .01  # fraction of objects to browse. Set to 1 to look at them all
    ncand = len(data['tile_id'])
    nfrac = int(fracbrowse * ncand)
    flag = np.zeros(nfrac)
    ids = np.random.permutation(np.arange(ncand))[0:nfrac]
    
    print '0: bad; 1: ok; 2: excellent!'
    for iw in range(ncand):
      objid = data['object_id']
      tileid = data['tile_id']
      url = 'http://upad.cefca.es/catalogues/jplus-test-2/object_query.html?image=%d&number=%d' % (tileid[ids[iw]],objid[ids[iw]])
      wb.open(url,new=0)
      flag[iw] = input('Flag object?: ')

    
