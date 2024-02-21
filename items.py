import math
import numpy as np
import spekpy as sp

class tube: 
  dTD = 100 #cm: The distance at which the fluence is calculated for the spectrum 
  coneangle = 20 #deg
  filters=[]
  def __init__(self, kV, mA, filters=[], dkV=0.5):
    self.s = sp.Spek(kvp=kV, mas = mA, th=20, dk=dkV, z=self.dTD, targ='W', char=False) # Generate a spectrum without characteristic Xray freq, and a Tungsten target at a 20deg angle. Fluence is per second (mas=mA)
    self.filters.append(('Be',0.8))
    for i in filters:
      self.filters.append(i)
    self.create_spectrum() ;
  def create_spectrum(self):
    print("Filters: ", self.filters)
    self.s.multi_filter(self.filters) # Filter by 0.8 mm of Be
    [self.freq, self.fluence] = self.s.get_spectrum()
    self.dk = self.freq[1]-self.freq[0]
    #minfreq = 40kV
    minfreq = np.argmax(self.freq>40)
    self.freq=self.freq[minfreq:] 
    self.fluence=self.fluence[minfreq:]
  def in_cone(self, px):
    if np.arccos(px[2]/np.linalg.norm(px)) < self.coneangle*np.pi/180:
      return True
    else:
      return False 
  

class detector:
  tube_distance = 1000 #mm
  offset_horiz = 0
  offset_vert = 0
  px_pitch = 1.27 #mm
  px_area = 0.0127*0.0127 #cm^2
  horiz_px = 200
  vert_px = 100
  def extent(self) : 
    return [self.offset_horiz-self.horiz_px/2*self.px_pitch, self.offset_horiz+self.horiz_px/2*self.px_pitch,
            self.offset_vert -self.vert_px /2*self.px_pitch, self.offset_vert +self.vert_px /2*self.px_pitch]

#======================= OBJECTS ================================================
class obj_sphere:
  def __init__ (self, center, radius, material, density, kv): #all distances in mm, density in g.cm^3
    self.material = material ; 
    self.density = density
    table = np.genfromtxt(material+".dat", skip_header=1)
    self.absorption = np.interp(kv, table[:,0], table[:,1]) * density 
    self.center = center ; 
    self.radius = radius 
  def ray_intersect_length(self, pxloc):
    a=(pxloc[0]**2+pxloc[1]**2+pxloc[2]**2)
    b=-2*(pxloc[0]*self.center[0]+pxloc[1]*self.center[1]+pxloc[2]*self.center[2])
    c=self.center[0]**2+self.center[1]**2+self.center[2]**2-self.radius**2
    delta = b**2-4*a*c
    if delta<=0: return 0
    alpha1 = (-b - np.sqrt(delta))/(2*a)
    alpha2 = (-b + np.sqrt(delta))/(2*a)
    return abs((alpha2-alpha1)*np.sqrt(pxloc[0]**2+pxloc[1]**2+pxloc[2]**2))*0.1 #for cm of intersection 

class obj_box:  
  def __init__ (self, center, size, material, density, kv):
    self.material = material ; 
    self.density = density
    table = np.genfromtxt(material+".dat", skip_header=1)
    self.absorption = np.interp(kv, table[:,0], table[:,1]) * density
    
    self.orig1 = np.array([center[0]-size[0]/2,center[1]-size[1]/2,center[2]-size[2]/2])
    self.orig2 = np.array([center[0]+size[0]/2,center[1]+size[1]/2,center[2]+size[2]/2])
    self.size = size
    self.planes = [subobj_rect(self.orig1, self.orig1+[size[0],0,0], self.orig1+[0,size[1],0]), 
                   subobj_rect(self.orig1, self.orig1+[size[0],0,0], self.orig1+[0,0,size[2]]),
                   subobj_rect(self.orig1, self.orig1+[0,size[1],0], self.orig1+[0,0,size[2]]),
                   subobj_rect(self.orig2, self.orig2-[size[0],0,0], self.orig2-[0,size[1],0]),
                   subobj_rect(self.orig2, self.orig2-[size[0],0,0], self.orig2-[0,0,size[2]]),
                   subobj_rect(self.orig2, self.orig2-[0,size[1],0], self.orig2-[0,0,size[2]]),]
  def ray_intersect_length(self,pxloc) :
    inter_pts = []
    for i in self.planes:
      inter = i.ray_intersect_pt(pxloc) ; 
      if len(inter)>0:
        inter_pts.append(inter) ; 
    if len(inter_pts)==2:
      return np.linalg.norm(inter_pts[0]-inter_pts[1])*0.1 #*0.1 to be in cm
    elif len(inter_pts) != 0:
      print(f"What the? Can't have {len(inter_pts)} intersection point with a box?")
    return 0       

class obj_cylinder:
  def __init__ (self, center, radius, height, material, density, kv):
    self.material = material ; 
    self.density = density
    table = np.genfromtxt(material+".dat", skip_header=1)
    self.absorption = np.interp(kv, table[:,0], table[:,1]) * density
    
    self.origin = np.array(center)
    self.radius = radius
    self.height = height
    self.ends = [subobj_disk(self.origin - np.array([0, self.height/2, 0]), [0,-1,0], self.radius), 
                 subobj_disk(self.origin + np.array([0, self.height/2, 0]), [0,1,0], self.radius)]    
  def ray_intersect_length (self,pxloc):
    inter_pts = []    
    a=(pxloc[0]**2+pxloc[2]**2)
    b=-2*(pxloc[0]*self.origin[0]+pxloc[2]*self.origin[2])
    c=self.origin[0]**2+self.origin[2]**2-self.radius**2
    delta = b**2-4*a*c
    if delta<=0: return 0
    alpha1 = (-b - np.sqrt(delta))/(2*a)
    alpha2 = (-b + np.sqrt(delta))/(2*a)
    
    if alpha1*pxloc[1] > self.origin[1]-self.height/2 and alpha1*pxloc[1] < self.origin[1]+self.height/2:
      inter_pts.append(alpha1 * pxloc) 
    if alpha2*pxloc[1] > self.origin[1]-self.height/2 and alpha2*pxloc[1] < self.origin[1]+self.height/2:
      inter_pts.append(alpha2 * pxloc)   
    
    if len(inter_pts) == 2:
      return np.linalg.norm(inter_pts[0]-inter_pts[1])*0.1
    
    for i in self.ends:
      inter = i.ray_intersect_pt(pxloc) ; 
      if len(inter)>0:
        inter_pts.append(inter) ; 
    if len(inter_pts)==2:
      return np.linalg.norm(inter_pts[0]-inter_pts[1])*0.1 #*0.1 to be in cm
    elif len(inter_pts) != 0:
      print(f"What the? Can't have {len(inter_pts)} intersection point with a cylinder?")
    return 0       
    
class subobj_disk:
  def __init__ (self, center, normal, radius):
    self.center = center 
    self.normal = normal
    self.radius = radius
  def ray_intersect_pt(self, px):
    if np.dot(px, self.normal) == 0:
      alpha = np.dot(self.center, px)
    else:
      alpha = np.dot(self.center, self.normal) / np.dot(px, self.normal)
    pt = alpha * px 
    t = pt - self.center
    if np.linalg.norm(t) < self.radius:
      return pt
    else:
      return []
    
class subobj_rect:
  def __init__(self, origin, pt1, pt2):
    self.v1 = pt1-origin
    self.v2 = pt2-origin
    self.origin = origin
    self.normal = np.cross(self.v1, self.v2)
    self.normal /= np.linalg.norm(self.normal)
    if abs(np.dot(self.v1, self.v2))>1e-7:
      print("The subobj_rect origin and point provided do not form orthogonal vectors")
  def ray_intersect_pt (self, px):
    if np.dot(px, self.normal) == 0:
      alpha = np.dot(self.origin, px)
    else:
      alpha = np.dot(self.origin, self.normal) / np.dot(px, self.normal)
    pt = alpha * px 
    t = pt - self.origin
    check1 = np.dot(self.v1, t) / np.dot(self.v1, self.v1)
    check2 = np.dot(self.v2, t) / np.dot(self.v2, self.v2)
    if (check1 >= 0 and check2 >= 0 and check1 <= 1 and check2 <= 1):
      return pt
    else :
      return [] 
  
#=======================================================================
class material:
  def __init__(self, name, kv):
    self.name = name
    self.table = np.genfromtxt(name+".dat", skip_header=1)
    self.absorptions = np.interp(kv, self.table[:,0], self.table[:,1])
  
#=======================================================================
def image (S, D, objects):
  img = np.zeros((D.horiz_px, D.vert_px, len(S.freq)))
  for i in range(0,D.horiz_px):
    for j in range(0,D.vert_px):
      pxloc = np.array([D.offset_horiz + (i-D.horiz_px/2)*D.px_pitch, D.offset_vert + (j-D.vert_px/2)*D.px_pitch, D.tube_distance])
      if not S.in_cone(pxloc):
        img[i,j,:] = 0
        continue 
      squarelaw = S.dTD**2/(np.dot(pxloc,pxloc)/100) #cm units
      img[i,j,:] = S.fluence * squarelaw
      for o in objects:
        length = o.ray_intersect_length(pxloc)
        img[i,j,:] *= np.exp(-o.absorption*length) ; 

  img_integr = np.sum(img, axis=-1)* S.dk * D.px_area
  return [img, img_integr] ; 




