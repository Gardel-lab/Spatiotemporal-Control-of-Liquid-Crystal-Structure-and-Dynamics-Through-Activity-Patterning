#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 10:17:29 2020

@author: stevenredford
"""

import numpy as np
from matplotlib import pyplot as plt
from numpy import genfromtxt


npts = 0
max_skip = 7
min_skip = 2
to_skip = 2
sincos = 1
nbins = 16

outlier_flag = 0

if outlier_flag == 1:    
    ntrks = 4
    upper_lim = 0.035
else:
    ntrks = 5
    upper_lim = 0.03
    

for nskip in range(min_skip,max_skip,to_skip):
    
    npts = npts + 1
    
    thetaz = []
    
    root = '/Users/stevenredford/Dropbox/tracks/'
    
    the_bins = np.linspace(0,180,nbins)
    
    
    for ii in range(ntrks):
        if ii == 0:    
            tail = 'SR1145B_2nd.csv'
        elif ii == 1:
            tail = 'SR1146A_1st.csv'
        elif ii == 2:
            tail = 'SR1146C_track.csv'
        elif ii == 3:
            tail = 'SR1147A_defect_direction.csv'
        elif ii == 4:
            tail = 'deflected_trajectory_1.csv'
        
        path = root + tail
        scale_fac = 0.1625
        tracks = genfromtxt(path, delimiter = ',')
        
        """
        if nskip == 2:
            plt.plot(tracks[:,3],tracks[:,4])
            plt.show()
        """
        
        subsampled = tracks[0::nskip,:]
        shp = np.shape(subsampled)

        
        for tt in range(shp[0] - 2):
            u = subsampled[tt + 1,3:5] - subsampled[tt,3:5]
            v = subsampled[tt + 2,3:5] - subsampled[tt + 1,3:5]
            
            
            #print('u = ', np.linalg.norm(u))
            #print('v = ', np.linalg.norm(v))
            
            if sincos == 0:    
                theta = np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))
            else:
                theta = np.arcsin(np.cross(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))
            thetaz.append(theta)
        
    
    
    thetaz = np.asarray(thetaz)*180/np.pi
    thetaz = thetaz[~np.isnan(thetaz)]
    
    globals()['theta_directed_%s' % nskip] = thetaz


    """
    If you wish to compare more data, this will import data into 'undirected'
    only 'directed' and 'bulk' are considered in the manuscript
    """
    root = '/Users/stevenredford/Dropbox/tracks/'
    #nskip = 2
    #nbins = 15
    
    thetas = []
    
    
    tail2 = '151A_three_proc.csv'
        
    path2 = root + tail2
    scale_fac = 0.1083
    
    tracks2 = genfromtxt(path2, delimiter = ',')
    tracks = tracks2
    
    dims = np.shape(tracks)
    
    track_labels = np.unique(tracks[:,1])
    ntracks = np.size(track_labels)
    #print(ntracks)
    
    for ii in track_labels:
        single_track = tracks[np.where(tracks[:,1] == ii),:]
        single_track = np.squeeze(single_track)
        shp_trk = np.shape(single_track)
        
        subsampled = single_track[0::nskip,:]
        shp = np.shape(subsampled)
        #print(shp)
        #if shp[0] > 20:
        #    continue
    
        for tt in range(shp[0] - 2):
            u = subsampled[tt + 1,3:5] - subsampled[tt,3:5]
            v = subsampled[tt + 2,3:5] - subsampled[tt + 1,3:5]
            
            #u = subsampled[1,3:5] - subsampled[0,3:5]
            #v = subsampled[tt + 2,3:5] - subsampled[tt + 1,3:5]
            
            if sincos == 0:    
                theta = np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))
            else:
                theta = np.arcsin(np.cross(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))
            thetas.append(theta)
        
    
    
    thetas = np.asarray(thetas)*180/np.pi
    thetas = thetas[~np.isnan(thetas)]
    
    globals()['theta_undirected_%s' % nskip] = thetas
    """
    """
    
    
    
    root = '/Users/stevenredford/Dropbox/tracks/'
    #nskip = 2
    #nbins = 15
    
    thetaaa = []
    
    tail2 = 'SR1129b_first_GS_proc.csv'
        
    path2 = root + tail2
    scale_fac = 0.1625
    
    tracks2 = genfromtxt(path2, delimiter = ',')
    tracks = tracks2
    #tracks = tracks2[np.where(tracks2[:,2] > 22),:]
    #tracks = np.squeeze(tracks)
    
    dims = np.shape(tracks)
    
    track_labels = np.unique(tracks[:,1])
    ntracks = np.size(track_labels)
    #print(ntracks)
    
    for ii in track_labels:
        single_track = tracks[np.where(tracks[:,1] == ii),:]
        single_track = np.squeeze(single_track)
        shp_trk = np.shape(single_track)
        """
        try:
            assert (shp_trk[1] == 8)
        except:
            continue
        """
        #print(shp_trk)
        subsampled = single_track[0::nskip,:]
        shp = np.shape(subsampled)
        #print(shp)
        #if shp[0] > 20:
        #    continue
    
        for tt in range(shp[0] - 2):
            u = subsampled[tt + 1,3:5] - subsampled[tt,3:5]
            v = subsampled[tt + 2,3:5] - subsampled[tt + 1,3:5]
            
            #u = subsampled[1,3:5] - subsampled[0,3:5]
            #v = subsampled[tt + 2,3:5] - subsampled[tt + 1,3:5]
            
            if sincos == 0:    
                theta = np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))
            else:
                theta = np.arcsin(np.cross(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))
            thetaaa.append(theta)
        
    
    
    thetaaa = np.asarray(thetaaa)*180/np.pi
    thetaaa = thetaaa[~np.isnan(thetaaa)]
    
    globals()['theta_bulk_%s' % nskip] = thetaaa
    

if sincos == 0:    
    the_bins = np.linspace(0,180,nbins)
else:    
    the_bins = np.linspace(-90,90,nbins)



#fig, ax = plt.subplots()
#ax.tick_params(axis="y",direction="in")
#ax.tick_params(axis="x",direction="in")
for ll in range(min_skip,max_skip,to_skip):
    temp = globals()['theta_directed_%d' % ll]
    #c = 0.25*(ll - 1)
    c = 0
    #a = 1 - ((ll - 1)*0.1)
    a = 1
    fig, ax = plt.subplots()
    ax.tick_params(axis="y",direction="in")
    #ax.tick_params(axis="x",direction="in")
    if sincos == 0:    
        plt.hist(temp, bins = the_bins, density = True, color = plt.cm.gray(c), edgecolor = plt.cm.gray(0.5), label = 'delta = %d' %ll, alpha = a)
    else:
        plt.hist(temp, bins = the_bins, density = True, color = plt.cm.gray(c), edgecolor = plt.cm.gray(0.5), label = 'delta = %d' %ll, alpha = a)
    plt.ylim([0, upper_lim])
    plt.legend() 
    plt.axvline(0, c = 'gray', ls = '--')
    plt.title('directed')
    pth = '/Users/stevenredford/Dropbox/figure_drafts_and_panels/turn_angle_asthetic_directed_%d' %ll 
    pth = pth + '.eps'
    #plt.savefig(pth)
    plt.show()    

"""
fig, ax = plt.subplots()
ax.tick_params(axis="y",direction="in")
ax.tick_params(axis="x",direction="in")
for ll in range(min_skip,max_skip,to_skip):
    temp = globals()['theta_undirected_%d' % ll]
    c = 0.25*(ll - 1)
    a = 1 - ((ll - 1)*0.1)
    if sincos == 0:    
        plt.hist(temp, bins = the_bins, density = True, color = plt.cm.gray(c), edgecolor = plt.cm.gray(0), label = 'delta = %d' %ll, alpha = a)
    else:
        plt.hist(temp, bins = the_bins, density = True, color = plt.cm.gray(c), edgecolor = plt.cm.gray(0), label = 'delta = %d' %ll, alpha = a)
plt.ylim([0, upper_lim])
plt.legend() 
plt.title('undirected')
#plt.savefig('/Users/stevenredford/Dropbox/figure_drafts_and_panels/turn_angle_undirected_050520.eps')
plt.show()    


#fig, ax = plt.subplots()
#ax.tick_params(axis="y",direction="in")
#ax.tick_params(axis="x",direction="in")
"""
for ll in range(min_skip,max_skip,to_skip):
    temp = globals()['theta_bulk_%d' % ll]
    #c = 0.25*(ll - 1)
    c = 0
    #a = 1 - ((ll - 1)*0.1)
    a = 1
    fig, ax = plt.subplots()
    ax.tick_params(axis="y",direction="in")
    #ax.tick_params(axis="x",direction="in")
    if sincos == 0:    
        plt.hist(temp, bins = the_bins, density = True, color = plt.cm.gray(c), edgecolor = plt.cm.gray(0.5), label = 'delta = %d' %ll, alpha = a)
    else:
        plt.hist(temp, bins = the_bins, density = True, color = plt.cm.gray(c), edgecolor = plt.cm.gray(0.5), label = 'delta = %d' %ll, alpha = a)
    plt.ylim([0, upper_lim])
    plt.legend() 
    plt.axvline(0, c = 'gray', ls = '--')
    plt.title('bulk')
    pth = '/Users/stevenredford/Dropbox/figure_drafts_and_panels/turn_angle_asthetic_unirected_%d' %ll
    pth = pth + '.eps'
    #plt.savefig(pth)
    plt.show()    
    


    
    
    
    

