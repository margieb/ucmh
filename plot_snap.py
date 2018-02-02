# open file
with open(filename,'rb') as ifile:
  delta = np.fromfile(ifile,count=-1,dtype=np.float32)
  delta.shape = (ncell,ncell,ncell)

# project over one of the dimensions
delta = np.average(delta,axis=2)

# set x and y scale
scale = (np.arange(ncell)+0.5)*boxsize*1./gridsize

# make 2d color plot with color bar
fig, ax = plt.subplots(figsize=(12.9, 10))
pcm = ax.pcolormesh(scale,scale,delta.T,cmap='cubehelix')
cbar = fig.colorbar(pcm, ax=ax)