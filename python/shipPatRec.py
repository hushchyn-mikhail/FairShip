__author__ = 'Mikhail Hushchyn'

import numpy as np

# Globals
ReconstructibleMCTracks = []
theTracks = []


def initialize(fgeo):
    pass


def execute(SmearedHits):
    """
    Main function of track pattern recognition.
    
    Parameters:
    -----------
    SmearedHits : list
        Smeared hits. SmearedHits = [{'digiHit': key, 
                                      'xtop': xtop, 'ytop': ytop, 'z': ztop, 
                                      'xbot': xbot, 'ybot': ybot, 
                                      'dist': dist2wire, 'detID': detID}, {...}, ...]
    """

    # withDist2Wire = False
    
    fittedtrackids = []
    track_hits = {}
    if len(SmearedHits) > 500:
        print "Too large hits in the event!"
        return track_hits
    
    # Separate hits
    SmearedHits_12y = []
    SmearedHits_12stereo = []
    SmearedHits_34y = []
    SmearedHits_34stereo = []
    for i_hit in range(len(SmearedHits)):
        ahit = SmearedHits[i_hit]
        detID = ahit['detID']
        statnb, vnb, pnb, lnb, snb = decodeDetectorID(detID)
        is_y12 = ((statnb == 1) + (statnb == 2)) * ((vnb == 0) + (vnb == 3))
        is_stereo12 = ((statnb == 1) + (statnb == 2)) * ((vnb == 1) + (vnb == 2))
        is_y34 = ((statnb == 3) + (statnb == 4)) * ((vnb == 0) + (vnb == 3))
        is_stereo34 = ((statnb == 3) + (statnb == 4)) * ((vnb == 1) + (vnb == 2))
        if is_y12:
            SmearedHits_12y.append(ahit)
        if is_stereo12:
            SmearedHits_12stereo.append(ahit)
        if is_y34:
            SmearedHits_34y.append(ahit)
        if is_stereo34:
            SmearedHits_34stereo.append(ahit)

    min_hits = 3
            
    #### PatRec in 12y
    short_tracks_y12 = pat_rec_view(SmearedHits_12y, min_hits)

    #### PatRec in 34y
    short_tracks_y34 = pat_rec_view(SmearedHits_34y, min_hits)

        
    # Extrapolation to the center of the magnet
    # ShipGeo is defined in macro/MufluxReco.py
    z_center = ShipGeo.Bfield.z

    # Combine track y12 and 34
    i_track_y12 = []
    i_track_y34 = []
    deltas_y = []
    for i_y12 in range(len(short_tracks_y12)):
        atrack_y12 = short_tracks_y12[i_y12]
        y_center_y12 = atrack_y12['k'] * z_center + atrack_y12['b']
        for i_y34 in range(len(short_tracks_y34)):
            atrack_y34 = short_tracks_y34[i_y34]
            y_center_y34 = atrack_y34['k'] * z_center + atrack_y34['b']
            i_track_y12.append(i_y12)
            i_track_y34.append(i_y34)
            deltas_y.append(abs(y_center_y12 - y_center_y34))

    max_dy = 50
    used_y12 = []
    used_y34 = []
    for i in np.argsort(deltas_y):
        dy = deltas_y[i]
        i_y12 = i_track_y12[i]
        i_y34 = i_track_y34[i]
        if dy < max_dy:
            if i_y12 not in used_y12:
                if i_y34 not in used_y34:
                    atrack = short_tracks_y34[i_y34]
                    atrack['i_track_y'] = i_y12
                    used_y12.append(i_y12)
                    used_y34.append(i_y34)

    ### PatRec in stereo12
    short_tracks_stereo12 = pat_rec_stereo_views(SmearedHits_12stereo, short_tracks_y12, min_hits)

    ### PatRec in stereo34
    short_tracks_stereo34 = pat_rec_stereo_views(SmearedHits_34stereo, short_tracks_y34, min_hits)

    # Prepare output of PatRec
    track_hits = {}
    for i_track in range(len(short_tracks_y12)):

        hits_y12 = []
        hits_stereo12 = []
        hits_y34 = []
        hits_stereo34 = []
        
        for ahit in short_tracks_y12[i_track]['hits']:
            hits_y12.append(ahit)
            
        for atrack in short_tracks_stereo12:
            if atrack['i_track_y'] == i_track:
                for ahit in atrack['hits']:
                    hits_stereo12.append(ahit)
                break

        for i_track_34 in range(len(short_tracks_y34)):
            atrack = short_tracks_y34[i_track_34]
            if atrack.has_key('i_track_y'):
                if atrack['i_track_y'] == i_track:
                    for ahit in atrack['hits']:
                        hits_y34.append(ahit)
                    break

        for atrack in short_tracks_stereo34:
            if atrack['i_track_y'] == i_track_34:
                for ahit in atrack['hits']:
                    hits_stereo34.append(ahit)
                break

        if len(hits_y12) >= min_hits and len(hits_stereo12) >= min_hits and len(hits_y34) >= min_hits and len(hits_stereo34) >= min_hits:
            atrack = {'y12': hits_y12, 'stereo12': hits_stereo12, 'y34': hits_y34, 'stereo34': hits_stereo34}
            track_hits[i_track] = atrack



    return track_hits


def finalize():
    pass


def decodeDetectorID(detID):
    """
    Decodes detector ID.

    Parameters
    ----------
    detID : int or array-like
        Detector ID values.

    Returns
    -------
    statnb : int or array-like
        Station numbers.
    vnb : int or array-like
        View numbers.
    pnb : int or array-like
        Plane numbers.
    lnb : int or array-like
        Layer numbers.
    snb : int or array-like
        Straw tube numbers.
    """

    statnb = detID // 10000000
    vnb = (detID - statnb * 10000000) // 1000000
    pnb = (detID - statnb * 10000000 - vnb * 1000000) // 100000
    lnb = (detID - statnb * 10000000 - vnb * 1000000 - pnb * 100000) // 10000
    snb = detID - statnb * 10000000 - vnb * 1000000 - pnb * 100000 - lnb * 10000 - 2000

    return statnb, vnb, pnb, lnb, snb


def hit_in_bin(x, y, k_bin, b_bin, k_size, b_size):
    """
        Counts hits in a bin of track parameter space (b, k).

        Parameters
        ---------
        x : array-like
            Array of x coordinates of hits.
        y : array-like
            Array of x coordinates of hits.
        k_bin : float
            Track parameter: y = k_bin * x + b_bin
        b_bin : float
            Track parameter: y = k_bin * x + b_bin

        Return
        ------
        track_inds : array-like
            Hit indexes of a track: [ind1, ind2, ...]
        """


    b_left = y - (k_bin - 0.5 * k_size) * x
    b_right = y - (k_bin + 0.5 * k_size) * x

    sel = (b_left >= b_bin - 0.5 * b_size) * (b_right <= b_bin + 0.5 * b_size) + \
    (b_left <= b_bin + 0.5 * b_size) * (b_right >= b_bin - 0.5 * b_size)

    return sel

def hit_in_window(x, y, k_bin, b_bin, window_width=1.):
    """
        Counts hits in a bin of track parameter space (b, k).

        Parameters
        ---------
        x : array-like
            Array of x coordinates of hits.
        y : array-like
            Array of x coordinates of hits.
        k_bin : float
            Track parameter: y = k_bin * x + b_bin
        b_bin : float
            Track parameter: y = k_bin * x + b_bin

        Return
        ------
        track_inds : array-like
            Hit indexes of a track: [ind1, ind2, ...]
        """


    y_approx = k_bin * x + b_bin

    flag = False
    if np.abs(y_approx - y) <= window_width:
        flag = True

    return flag


def get_zy_projection(z, xtop, ytop, xbot, ybot, k_y, b_y):
    
    x = k_y * z + b_y
    
    k = (ytop - ybot) / (xtop - xbot + 10**-6)
    b = ytop - k * xtop
    y = k * x + b
    
    return y


def pat_rec_view(SmearedHits, min_hits):

    long_tracks = []

    # Take 2 hits as a track seed
    for ahit1 in SmearedHits:
        for ahit2 in SmearedHits:

            if ahit1['z'] >= ahit2['z']:
                continue
            if ahit1['detID'] == ahit2['detID']:
                continue

            # +- dist2wire
            for sign1 in [-1, 1]:
                for sign2 in [-1, 1]:

                    y1 = ahit1['ytop'] + sign1 * ahit1['dist']
                    y2 = ahit2['ytop'] + sign2 * ahit2['dist']
                    z1 = ahit1['z']
                    z2 = ahit2['z']
                    layer1 = ahit1['detID'] // 10000
                    layer2 = ahit2['detID'] // 10000

                    k_bin = 1. * (y2 - y1) / (z2 - z1)
                    b_bin = y1 - k_bin * z1

                    if abs(k_bin) > 1:
                        continue

                    atrack = {}
                    atrack['hits'] = [ahit1, ahit2]
                    atrack['y'] = [y1, y2]
                    atrack['z'] = [z1, z2]
                    atrack['layer'] = [layer1, layer2]

                    # Add new hits to the seed
                    for ahit3 in SmearedHits:

                        if ahit3['detID'] == ahit1['detID'] or ahit3['detID'] == ahit2['detID']:
                            continue

                        dist3 = ahit3['dist']
                        y3 = ahit3['ytop']
                        z3 = ahit3['z']
                        layer3 = ahit3['detID'] // 10000

                        if layer3 in atrack['layer']:
                            continue

                        in_bin_p = hit_in_window(z3, y3+dist3, k_bin, b_bin, window_width=.2)
                        err_p = abs(k_bin * z3 + b_bin - (y3+dist3))
                        in_bin_m = hit_in_window(z3, y3-dist3, k_bin, b_bin, window_width=.2)
                        err_m = abs(k_bin * z3 + b_bin - (y3-dist3))
                        if in_bin_p or in_bin_m:
                            atrack['hits'].append(ahit3)
                            atrack['z'].append(z3)
                            atrack['layer'].append(layer3)
                            if err_m < err_p:
                                atrack['y'].append(y3-dist3)
                            else:
                                atrack['y'].append(y3+dist3)

                    if len(atrack['hits']) >= min_hits:
                        long_tracks.append(atrack)

    # Remove clones
    used_hits = []
    short_tracks = []
    n_hits = [len(atrack['hits']) for atrack in long_tracks]

    for i_track in np.argsort(n_hits)[::-1]:

        atrack = long_tracks[i_track]
        new_track = {}
        new_track['hits'] = []
        new_track['y'] = []
        new_track['z'] = []
        new_track['layer'] = []

        for i_hit in range(len(atrack['hits'])):
            ahit = atrack['hits'][i_hit]
            if ahit['digiHit'] not in used_hits:
                new_track['hits'].append(ahit)
                new_track['y'].append(atrack['y'][i_hit])
                new_track['z'].append(atrack['z'][i_hit])

        if len(new_track['hits']) >= min_hits:
            short_tracks.append(new_track)
            for ahit in new_track['hits']:
                used_hits.append(ahit['digiHit'])

    # Fit tracks
    for atrack in short_tracks:
        [atrack['k'], atrack['b']] = np.polyfit(atrack['z'], atrack['y'], deg=1)

    return short_tracks


def pat_rec_stereo_views(SmearedHits_stereo, short_tracks_y, min_hits):

    ### PatRec in stereo12
    long_tracks_stereo = []
    short_tracks_stereo = []
    used_hits = []
    deg = np.deg2rad(ShipGeo.strawtubes.ViewAngle)

    for i_track_y in range(len(short_tracks_y)):

        atrack_y = short_tracks_y[i_track_y]
        k_y = atrack_y['k']
        b_y = atrack_y['b']

        # Get hit zx projections
        for ahit in SmearedHits_stereo:
            x_center  = get_zy_projection(ahit['z'], ahit['ytop'], ahit['xtop'],ahit['ybot'], ahit['xbot'], k_y, b_y)
            ahit['zx_projection'] = x_center

        temp_tracks_stereo = []

        for ahit1 in SmearedHits_stereo:
            for ahit2 in SmearedHits_stereo:

                if ahit1['z'] >= ahit2['z']:
                    continue
                if ahit1['detID'] == ahit2['detID']:
                    continue
                if ahit1['digiHit'] in used_hits:
                    continue
                if ahit2['digiHit'] in used_hits:
                    continue

                x1_center  = ahit1['zx_projection']
                x2_center  = ahit2['zx_projection']

                if abs(x1_center ) > 300 or abs(x2_center ) > 300:
                    continue

                for sign1 in [-1, 1]:
                    for sign2 in [-1, 1]:

                        x1 = x1_center + sign1 * ahit1['dist'] / np.sin(deg)
                        x2 = x2_center + sign2 * ahit2['dist'] / np.sin(deg)
                        z1 = ahit1['z']
                        z2 = ahit2['z']
                        layer1 = ahit1['detID'] // 10000
                        layer2 = ahit2['detID'] // 10000

                        k_bin = 1. * (x2 - x1) / (z2 - z1)
                        b_bin = x1 - k_bin * z1

                        atrack = {}
                        atrack['hits'] = [ahit1, ahit2]
                        atrack['x'] = [x1, x2]
                        atrack['z'] = [z1, z2]
                        atrack['layer'] = [layer1, layer2]
                        atrack['i_track_y'] = i_track_y

                        for ahit3 in SmearedHits_stereo:

                            if ahit3['digiHit'] == ahit1['digiHit'] or ahit3['digiHit'] == ahit2['digiHit']:
                                continue
                            if ahit3['digiHit'] in used_hits:
                                continue

                            dist3 = ahit3['dist'] / np.sin(deg)
                            x3_center = ahit3['zx_projection']
                            z3 = ahit3['z']
                            layer3 = ahit3['detID'] // 10000

                            if abs(x3_center) > 300:
                                continue

                            if layer3 in atrack['layer']:
                                continue

                            in_bin_p = hit_in_window(z3, x3_center+dist3, k_bin, b_bin, window_width=.5)
                            err_p = abs(k_bin * z3 + b_bin - (x3_center+dist3))
                            in_bin_m = hit_in_window(z3, x3_center-dist3, k_bin, b_bin, window_width=.5)
                            err_m = abs(k_bin * z3 + b_bin - (x3_center-dist3))
                            if in_bin_p or in_bin_m:
                                atrack['hits'].append(ahit3)
                                atrack['z'].append(z3)
                                atrack['layer'].append(layer3)
                                if err_m < err_p:
                                    atrack['x'].append(x3_center-dist3)
                                else:
                                    atrack['x'].append(x3_center+dist3)

                        if len(atrack['hits']) >= min_hits:
                            temp_tracks_stereo.append(atrack)
                            long_tracks_stereo.append(atrack)


        # Remove clones
        max_track = None
        max_n_hits = -999

        for atrack in temp_tracks_stereo:
            if len(atrack['hits']) > max_n_hits:
                max_track = atrack
                max_n_hits = len(atrack['hits'])

        if max_track is not None:
            short_tracks_stereo.append(max_track)
            for ahit in max_track['hits']:
                used_hits.append(ahit['digiHit'])

    return short_tracks_stereo