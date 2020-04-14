from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import numpy as np
import sys
cubelist=open('fields.txt','r')
used=0
wave_step=0.5
wave_start=4600.0/(1+1.5)
wave_end=8000
    
dest_wave=np.arange(wave_start, wave_end, wave_step)
dest_count=[]
dest_n=[]
for i in range(len(dest_wave)):
        dest_count.append(0)
        dest_n.append(0)
        
print wave_start, wave_end, wave_step, len(dest_wave)
    
for cubeid in cubelist:
    if "#" in cubeid: continue
    if "X" in cubeid: break
    qso=cubeid.split()[0]
    
    #orginal fluxcube
    #qso in fields.txt holds all identifiers such as 'J0014m0028'
    cube=Cube('%s/cube.fits' % qso)
    
    #the (current) final full catalogs
    #http://tucana.astro.physik.uni-potsdam.de/~mwendt/2020/mf/catalogs_final/
    catalog=open('../mfcatalogs/catalog_%s.txt' % qso,'r')
    tags=catalog.readline()
    average_ptr=tags.split(',').index(' average')
    print 'pointer to average:',average_ptr
    print 'qso', qso
    
    plt.figure(figsize=(16,9))
    count=0
    

    for line in catalog:
        if 'id' in line: continue
        #id, px, py, z, b, quality, template, used, gaussx, gaussy, gaussrot, smoothratio, oiiratio, qso_ra, qso_dec, qso_x, qso_y, average, sigma, median, martin, johannes, nicolas
        px=float(line.split(',')[1])
        py=float(line.split(',')[2])
        pz=float(line.split(',')[3])
        q =int(line.split(',')[5])
        average= float(line.split(',')[average_ptr])
        
        #selection criteria here
        if (average<1.5):
        #if (average>1.5 or average<0.8):
            continue
       
        ds=3 #spatial size is +/- ds
        px=int(px)
        py=int(py)
        print px,py,pz,q
        spec=cube[:,int(py)-ds:int(py)+ds,int(px)-ds:int(px)+ds].mean(axis=(1,2))
        
        data0=spec.data
        data1=np.nan_to_num(data0)
        used+=1
        
        #use debug spectrum instead
        #for i in range(len(data1)):
        #    spec.data[i]=0
        #    if i> count*100: 
        #        spec.data[i]=1
        #pz=0#count 
            
        waven=np.arange(spec.wave.get_crval(), spec.wave.get_crval()+spec.wave.get_step()*spec.wave.shape,spec.wave.get_step())
        
        #use debug spectrum instead
        #waven=np.arange(0,0+len(data1)*spec.wave.get_step(),spec.wave.get_step())
            
        #translate to restframe
        for i in range(len(waven)):
            waven[i] /= (pz+1)
            
        #populate destination spectrum
        
        #estimate starting point
        first = waven[0]
        
        if first < wave_start:
            kfrom=0
        else: #spec starts later than template
            kfrom = int ( (first - wave_start) / wave_step - 5)
        
        
        
        mfrom=0
        mto=len(waven)
        for k in range(kfrom, len(dest_wave)):
            #d destination
            #s start
            #e end
            #r range
            #c count
            dr = wave_step
            ds = wave_start +  wave_step*k - wave_step*0.5
            de = wave_start +  wave_step*k + wave_step*0.5
            #print ds
            hit=False
            for m in range(mfrom,mto):
                #a slight local approximation as m-1 and m+1 to m differ
                if (m < (len(waven)-1)):
                    dm = waven[m+1]-waven[m]
                ms = waven[m] - dm*0.5
                me = waven[m] + dm*0.5
                mr = me-ms
                
                #no overlap
                if hit==False:
                    if ms > de or me < ds: continue
                
                if hit==True: #we already filled the new pixel and are beyond it
                    if ms > de or me < ds:
                        #mto=min(len(waven),m+5)
                        break
                
                mfrom=max(0,m-5)
                mto=min(m+10,len(waven))
                hit=True
                mc = data1[m]
                
                #pixel fully included in destination range 
                if ms >= ds and me <= de:
                    if not mr < dr:
                        print 'error 1'
                    dest_count[k] +=  mc
                    dest_n[k]+= 1 #(mr/dr)
                    
                #full overlap    
                if ms < ds and me > de:
                    #add partially
                    frac = 1.0 #dr / mr 
                    dest_count[k] +=  (mc * frac)
                    dest_n[k]+=frac
                    
                #left overlap
                if ms <= ds and me <= de:
                    overlap = me - ds
                    frac = overlap / dr
                    dest_count[k] +=  (mc * frac)
                    dest_n[k]+=frac
                    
                #right overlap
                if ms >= ds and me >= de:
                    overlap = de - ms
                    frac = overlap / dr
                    dest_count[k] +=  (mc * frac)
                    dest_n[k]+=frac
        del spec
            
        #plt.set_xticks(np.arange(crval,crmax,200))
       
        #plt.step(dest_wave,dest_count,where='mid')
    del cube   
    
print used
plt.step(dest_wave,dest_count,where='mid')
plt.savefig('combined_all.png')

specfile=open('combined_all.asc' , 'w')

for i in range(len(dest_wave)):
    specfile.write('%.5f %.5f %.5f\n' % (dest_wave[i],dest_count[i],dest_n[i]))
    
