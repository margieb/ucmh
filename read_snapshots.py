with open("/nas02/home/m/n/mbruff/lustre/run1/output/snapshot_099", "rb") as f:
  # header
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be 256
  npart = np.fromfile(f,dtype=np.int32,count=6) # num particles of each type in file, DM particles are type 1
  massperparticle = np.fromfile(f,dtype=np.float64,count=6) # mass per particle in Msun/h
  scalefactor = np.fromfile(f,dtype=np.float64,count=1)
  redshift = np.fromfile(f,dtype=np.float64,count=1)
  f.seek(8,1) # skip 8 bytes
  ntotal = np.fromfile(f,dtype=np.int32,count=6) # num particles total (versus in this file)
  f.seek(4,1)
  numfiles = np.fromfile(f,dtype=np.int32,count=1)
  boxsize = np.fromfile(f,dtype=np.float64,count=1)
  omegaM = np.fromfile(f,dtype=np.float64,count=1)
  omegaL = np.fromfile(f,dtype=np.float64,count=1)
  hubbleparam = np.fromfile(f,dtype=np.float64,count=1)
  f.seek(96,1)
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be 256
  
  # positions
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be npart*12
  pos = np.fromfile(f.dtype=np.float32,count=npart[1]*3)
  pos.shape = (-1,3) # pos[i] gives x,y,z of i-th particle
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be npart*12
  
  # velocities
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be npart*12
  vel = np.fromfile(f.dtype=np.float32,count=npart[1]*3)
  vel.shape = (-1,3) # vel[i] gives vx,vy,vz of i-th particle
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be npart*12
  
  # id numbers (only useful if tracking a particle through multiple snapshots)
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be npart*4
  ids = np.fromfile(f.dtype=np.int32,count=npart[1])
  blockcount = np.fromfile(f,dtype=np.int32,count=1) # should be npart*4
