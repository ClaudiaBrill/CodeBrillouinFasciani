%%%%CODE FOR BRILLOUIN SPECTRA MEASUREMENTS%%
%%Used for Brillouin Microscopy measurements in the paper Fasciani et al., to appear in Nature Genetics

%%%Claudia Testi, PhD
%%Post-doc @ Istituto Italiano di Tecnologia
%%claudia.testi@iit.it

%%Fabrizio Gala, PhD
%%CrestOptics S.p.A.

%%%Please cite if used!

function  imFilt = smoothWithNeighbours2(im,nii,njj,treshold)

    imFilt = im ;
    
    nlist = 0 ;
    nn = zeros(nii(2), njj(2));
    for ii=nii(1):nii(2)
        for jj=njj(1):njj(2)
            if im(ii,jj)==treshold
                nlist = nlist+1 ;
                nn(ii,jj)=1 ; 
                list(nlist).ii  = ii;
                list(nlist).jj  = jj;
                list(nlist).noFit = true ;
            end
        end
    end
    
    counter = 0 ;
    cc  = 0 ;

    while(counter ~= nlist)

        for i=1:nlist
            
            if list(i).noFit == true
                ii = list(i).ii ;
                jj = list(i).jj ;

                cc = 0;
                mb = 0.0 ;
        
                if (ii+1)<= nii(2) && nn(ii+1,jj)==0  
                    cc = cc+1 ;
                    mb = mb + imFilt(ii+1,jj) ;
                end
                if (ii-1)>= nii(1) && nn(ii-1,jj)==0   
                    cc = cc+1 ;
                    mb = mb + imFilt(ii-1,jj) ;
                end
            
                if (jj+1)<= njj(2) && nn(ii,jj+1)==0    
                    cc = cc+1 ;
                    mb = mb + imFilt(ii,jj+1) ;
                end
                if (jj-1)>= njj(1) && nn(ii,jj-1)==0  
                    cc = cc+1 ;
                    mb = mb + imFilt(ii,jj-1) ;
                end
        
                if (ii+1)<= nii(2) && (jj-1) >= njj(1) && nn(ii+1,jj-1)==0 
                    cc = cc+1 ;
                    mb = mb + imFilt(ii+1,jj-1) ;
                end
                if (ii+1)<= nii(2) && (jj+1) <= njj(2) && nn(ii+1,jj+1)==0 
                    cc = cc+1 ;
                    mb = mb + imFilt(ii+1,jj+1) ;
                end
                if (ii-1)>= nii(1) && (jj-1) >= njj(1) && nn(ii-1,jj-1)==0 
                    cc = cc+1 ;
                    mb = mb + imFilt(ii-1,jj-1) ;
                end
                if (ii-1)>= nii(1) && (jj+1) <= njj(2) && nn(ii-1,jj+1)==0 
                    cc = cc+1 ;
                    mb = mb + imFilt(ii-1,jj+1) ;
                end
            
                if cc ~= 0
                    imFilt(ii,jj)= mb / cc;
                    nn(ii,jj)=0;
                    list(i).noFit = false;
            
                    counter = counter+1 ;
                end
            end
        end
    end
end

