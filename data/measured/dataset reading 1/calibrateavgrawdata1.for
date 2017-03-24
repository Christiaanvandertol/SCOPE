
	program calibrateavgrawdata

	implicit real(a-z)
      dimension sum(64),count(64),avg(64),channel(64)
      integer  i,j,line,nolines,day,daynr
	character header1,header2,header3,header4

        open (unit=3,file='rawplot1plot2data.txt',status='old')  
        open (unit=4,file='halfhourplot1and2.out',status='unknown') 
        open (unit=11,file='halfhourplot1and2from144.out',status='old') 
        open (unit=12,file='aws144-269.txt',status='old') 
        open (unit=13,file=
     .  'in-situandremoteGcombinedplusTC.txt',status='old') 
        open (unit=14,file='eddydatasonningmaize144-269.txt',status=
     . 'old') 
        open (unit=5,file='verivarsplot1.out',status='unknown')
        open (unit=6,file='verivarsplot2.out',status='unknown') 
	  open (unit=7,file='drivingvarsplot1.out',status='unknown')
        open (unit=8,file='drivingvarsplot2.out',status='unknown')    
        open (unit=9,file='maizebradley.dat',status='unknown')  
	   
        sigma = 5.67e-8

	  cal306=91.324
        cable306=2/0.385
        cal686=1000/13.4

	  cal307=81.766
        cable307=2/0.385
        cal687=1000/14.1



      read(11,'(a1)')header1
      read(12,'(a1)')header2
      read(13,'(a1)')header3
      read(14,'(a1)')header4

        goto 2  !this goes straight to the half-hourly file padded at the front from 144-162.
        nolines = 0

	do i =1,64 
          sum(i) = 0.0
        enddo

	do i =1,64 
          count(i) = 0.0
        enddo 

C-------106 days, from DOY 163-268
        do 1 line=1,15264
        nolines = nolines + 1
        read(3,*) (channel(j),j=1,64)
      

        do i=1,64
        if(channel(i).ne.-999) then
          sum(i) = sum(i) + channel(i)
          count(i) = count(i) + 1
        endif
        enddo
        
        if(nolines.eq.3) then
         do i=1,64
           if(count(i).gt.0.0) then
             avg(i) = sum(i)/count(i)
             sum(i) = 0.0
             count(i) = 0.0
           else
             avg(i) = -999
           endif
         enddo
        nolines = 0
          
c        write(4,'(64(f11.5,1x))') (avg(i),i=1,64)
        endif

1       continue

c        rewind(4)

c        do line =1,5088




2        do line =1,6000

	read(11,*) day_y,AN3,AN4_149cm,AN4_86cm,pp4022_1,pp4022_2,pp4022_3,
     .pp4022_4,pp4022_5,pp4022_6,irt4940,nrlit686,thetP15,thm95,thm105,
     .thm127,thm112,thm125,cmup307,cgup307,cmdn307,cgdn307,prt307,
     .RH_200,Tair_200,RH_400,Tair_400,windvane,RH_can,Tair_can,thm91,
     .thm87,thm129,thm33,thm41,thm15,thm10,thm13,thm6,thm58,cmup306,
     .cgup306,cmdn306,cgdn306,prt306,RH_c2,Tair_c2,thm35,pp4011_1,
     .pp4011_2,pp4011_3,pp4011_4,pp4011_5,pp4011_6,thetP24,irt4941,
     .nrlit687,nrlit680,thm5,thm12,thm4,thm19,thm16,SkyeP	

      write(*,*) Tair_can,Tair_400,Tair_c2
	read(12,*) DOYAWS,Sin_AWS,rain_AWS,u_aws,Ta_AWS,Tw_AWS,ea_AWS,
     .es_AWS,RH_AWS
  
      read(13,*) DOYG,G1_insitu,G2_insitu,NRlite680,NRlite682,NRlite698,
     .Gset40,Gset42,Gset43,Gset44,Gset45,Gset39,Gset46,
     .Gset38,G1_remoteavg, G1_remoterows,G1_remoteinter,Gset06,Gset08,
     .Gset07,Gset09,Gset20,Gset23,Gset24,Gset62,G2_remoteavg,
     .G2_remoterows,G2_remoteinter,sigmaG1,sigmaG1row,sigmaG1inter,
     .sigmaG2,sigmaG2row,sigmaG2inter,Tset06plot2R,Tset08plot2I,
     .Tset07plot2I,Tset09plot2R,Tset20plot2R,Tset23plot2I,Tset24plot2I,
     .Tset62plot2R,Tset40plot1I,Tset42plot1R,Tset43plot1I,Tset44plot1I,
     .Tset45plot1I,Tset39plot1R,Tset46plot1R,Tset38plot1R



	read(14,*) day_eddy,ZmHt,ZPD,windDir,u_eddy,T_eddy,q_eddy,
     .CO2_eddy,rhoCp,Monin,Ustr,ZoL,H_eddy,LE_eddy,Fco2			

c	write(*,*) day_y,DOYAWS,DOYG,day_eddy
     
C------Calibration of Profileprobes and Thetaprobes
C------Profileprobes plot 1
      if(pp4022_1.gt.0.0) then
	ppplot1_10=(0.65+13.6*(pp4022_1/1000)-29.9*(pp4022_1/1000)**2.0+
     .   72.5*(pp4022_1/1000)**3.0-1.6)/8.4
	else
	   ppplot1_10 = -999 
	endif
	
      if(pp4022_2.gt.0.0) then
	ppplot1_20=(0.65+13.6*(pp4022_2/1000)-29.9*(pp4022_2/1000)**2.0+
     .   72.5*(pp4022_2/1000)**3.0-1.6)/8.4
      else
	   ppplot1_20 = -999 
	endif
	
      if(pp4022_3.gt.0.0) then
	ppplot1_30=(0.65+13.6*(pp4022_3/1000)-29.9*(pp4022_3/1000)**2.0+
     .   72.5*(pp4022_3/1000)**3.0-1.6)/8.4
      else
	   ppplot1_30 = -999 
	endif
	
      if(pp4022_4.gt.0.0) then
	ppplot1_40=(0.65+(13.6*(pp4022_4/1000)-29.9*(pp4022_4/1000)**2.0+
     .   72.5*(pp4022_4/1000)**3.0)-1.6)/8.4
      else
 	   ppplot1_40 = -999 
	endif	
	
      if(pp4022_5.gt.0.0) then
	ppplot1_60=(0.65+13.6*(pp4022_5/1000)-29.9*(pp4022_5/1000)**2.0+
     .   72.5*(pp4022_5/1000)**3.0-1.6)/8.4
      else
	   ppplot1_60 = -999 
	endif
	
      if(pp4022_6.gt.0.0) then
	ppplot1_100=(0.65+13.6*(pp4022_6/1000)-29.9*
     .(pp4022_6/1000)**2.0+72.5*(pp4022_6/1000)**3.0-1.6)/8.4
      else
	   ppplot1_100 = -999 
	endif

C------Profileprobes plot 2
	
      if(pp4011_1.gt.0.0) then
	ppplot2_10=(0.65+13.6*(pp4011_1/1000)-29.9*(pp4011_1/1000)**2.0+
     .   72.5*(pp4011_1/1000)**3.0-1.6)/8.4
	else
	   ppplot2_10 = -999 
	endif
	
      if(pp4011_2.gt.0.0) then
	ppplot2_20=(0.65+13.6*(pp4011_2/1000)-29.9*(pp4011_2/1000)**2.0+
     .   72.5*(pp4011_2/1000)**3.0-1.6)/8.4
      else
	   ppplot2_20 = -999 
	endif
	
      if(pp4011_3.gt.0.0) then
	ppplot2_30=(0.65+13.6*(pp4011_3/1000)-29.9*(pp4011_3/1000)**2.0+
     .   72.5*(pp4011_3/1000)**3.0-1.6)/8.4
      else
	   ppplot2_30 = -999 
	endif
	
      if(pp4011_4.gt.0.0) then
	ppplot2_40=(0.65+(13.6*(pp4011_4/1000)-29.9*(pp4011_4/1000)**2.0+
     .   72.5*(pp4011_4/1000)**3.0)-1.6)/8.4
      else
 	   ppplot2_40 = -999 
	endif	
	
      if(pp4011_5.gt.0.0) then
	ppplot2_60=(0.65+13.6*(pp4011_5/1000)-29.9*(pp4011_5/1000)**2.0+
     .   72.5*(pp4011_5/1000)**3.0-1.6)/8.4
      else
	   ppplot2_60 = -999 
	endif
	
      if(pp4011_6.gt.0.0) then
	ppplot2_100=(0.65+13.6*(pp4011_6/1000)-29.9*
     .(pp4011_6/1000)**2.0+72.5*(pp4011_6/1000)**3.0-1.6)/8.4
      else
	   ppplot2_100 = -999 
	endif

C------Thetaprobes plot 1

      if(thetP15.gt.0.0) then	
        tpplot1_5 = (1.07+6.4*(thetP15/1000)-6.4*(thetP15/1000)**2.0+
     .  4.7*(thetP15/1000)**3.0-1.6)/8.4
	else
	  tpplot1_5 = -999
	endif 

C------Thetaprobes plot 2

      if(thetP24.gt.0.0) then     
        tpplot2_5 = (1.07+6.4*(thetP24/1000)-6.4*(thetP24/1000)**2.0+
     .  4.7*(thetP24/1000)**3.0-1.6)/8.4
      else
	  tpplot2_5 = -999
	endif
	
C-------Calculate radiation terms/surface temperatures

C-------Plot 1&2, driving variables

C-------Plot 1 NRlite, 6 m

        if(nrlit686.ne.-999) then
          rnnrliteplot1 = nrlit686*cal686
	  else
	    rnnrliteplot1 = -999
	  endif

C-------Plot 2 NRlite, 6 m

        if(nrlit687.ne.-999) then
          rnnrliteplot2 = nrlit687*cal687
	  else
	    rnnrliteplot2 = -999
	  endif

C-------Plot 1, 4 component, 6 m

        if(cmup306.ne.-999) then
          Sinplot1 = cal306*cmup306
	  else
          Sinplot1 = -999
	  endif
        if(cmdn306.ne.-999) then
          Soutplot1 = cal306*cmdn306
	  else
          Soutplot1 = -999
	  endif


C-------Plot 2, 4 component, 6 m

        if(cmup307.ne.-999) then
          Sinplot2 = cal307*cmup307
	  else
          Sinplot2 = -999
	  endif
        if(cmdn307.ne.-999) then
          Soutplot2 = cal307*cmdn307
	  else
          Soutplot2 = -999
	  endif
  
C-------to convert prt306 and prt307 to an outgoing longwave flux;


        if(prt306.ne.-999.and.cgup306.ne.-999.and.cgdn306.ne.-999) then
          Tprt306 = prt306 - cable306 + 273.15
          Linst_out306 = sigma*((Tprt306)**4.0)

          Linplot1 = cgup306*cal306+ Linst_out306 
          Loutplot1 =cgdn306*cal306+ Linst_out306 
	   else
	    Linplot1 = -999
          Loutplot1 = -999
	   endif

        if(prt307.ne.-999.and.cgup307.ne.-999.and.cgdn307.ne.-999) then
          Tprt307 = ((prt307-13.14)/1.08) - cable307 + 273.15
          Linst_out307 = sigma*((Tprt307)**4.0)

          Linplot2 = cgup307*cal307+ Linst_out307 
          Loutplot2 =cgdn307*cal307+ Linst_out307 
	   else
	    Linplot2 = -999
          Loutplot2 = -999
	   endif

         if(Sinplot1.ne.-999.and.Soutplot1.ne.-999.and.
     .   Linplot1.ne.-999.and.Loutplot1.ne.-999) then
           rn4compplot1 = Sinplot1-Soutplot1+Linplot1-Loutplot1
	   else
           rn4compplot1 = -999
	   endif

         if(Sinplot2.ne.-999.and.Soutplot2.ne.-999.and.
     .   Linplot2.ne.-999.and.Loutplot2.ne.-999) then
           rn4compplot2 = Sinplot2-Soutplot2+Linplot2-Loutplot2
	   else
           rn4compplot2 = -999
	   endif

         eps = 0.98
         eps = 0.95



	  if(Loutplot1.gt.0.0) then
          Ts4compmaizeplot1 = ((Loutplot1/sigma)**0.25)-273.15
          Ts4compmaizeplot1=(((Loutplot1-(1-eps)*Linplot1)/
     .    (eps*sigma))**0.25)-273.15
        else
          Ts4compmaizeplot1=-999
        endif

	  if(Loutplot2.gt.0.0) then
          Ts4compmaizeplot2 = ((Loutplot2/sigma)**0.25)-273.15
          Ts4compmaizeplot2=(((Loutplot2-(1-eps)*Linplot2)/
     .    (eps*sigma))**0.25)-273.15
        else
          Ts4compmaizeplot2=-999
        endif

	  if(irt4940.gt.0.0) then
         TsIRTplot1 = (irt4940/10.0)-0.59
	  else
	   TsIRTplot1 = -999
	  endif

	  if(irt4941.gt.0.0) then
         TsIRTplot2 = irt4941/10.0
	  else
	   TsIRTplot2 = -999
	  endif

C--------Longin from AWS mast

C---------Total sky emissivity according to Brutsaert

	    esky0_AWS= 1.24*((ea_AWS*10.0)/(Ta_AWS+273.15))**0.143
          Lin_AWS = esky0_AWS*sigma*(Ta_AWS+273.15)**4.0

	    


C-------Calculate actual vapour pressures
C-------Plot 1, driving variables

          if(RH_200.gt.0.0.and.RH_200.lt.100.and.Tair_200.gt.0.0) then
            ea_200 = es(Tair_200)*(RH_200/100.0)
	    else
            ea_200 = -999
          endif
          if(RH_400.gt.0.0.and.RH_400.lt.100.and.Tair_400.gt.0.0) then
            ea_400 = es(Tair_400)*(RH_400/100.0)
	    else
	      ea_400 = -999
	    endif

C-------Plot 1, verification variables

          if(RH_can.gt.0.0.and.RH_can.lt.100.and.Tair_can.gt.0.0) then
            ea_can = es(Tair_can)*(RH_can/100.0)
	    else
            ea_can = -999
          endif

C-------Plot 2, verification variables

          if(RH_c2.gt.0.0.and.RH_c2.lt.100.and.Tair_c2.gt.0.0) then
            ea_c2 = es(Tair_c2)*(RH_c2/100.0)
	    else
            ea_c2 = -999
          endif

c        write(7,'(22(f13.5,2x))') day_y,Tair_200,Tair_400,ea_200,ea_400,
c     .  ea_AWS*10,Sinplot1,Sin_AWS,Linplot1,Lin_AWS,u_eddy,u_aws

        write(7,'(22(f13.5,2x))') day_y,Tair_200,Tair_400,ea_200,ea_400,
     .  Sinplot1,Linplot1,u_aws,tpplot1_5,CO2_eddy

        write(8,'(22(f13.5,2x))') day_y,Tair_200,Tair_400,ea_200,ea_400,
     .  Sinplot2,Linplot2,u_aws,tpplot2_5,CO2_eddy

        write(9,'(22(f13.5,2x))') day_y,Tair_200,Tair_400,ea_200,ea_400,
     .  q_eddy,CO2_eddy,u_aws, Sinplot1,Linplot1,rn4compplot1,G1_insitu,
     .  TsIRTplot1,Ts4compmaizeplot1,H_eddy,LE_eddy,Fco2,Ustr,u_eddy,
     .  tpplot1_5



c- actual vapour pressure or vapour pressure deficit above the canopy
c- wind speed above the canopy
c- incoming radiation (spectrum, or at least integrated between ~0.4-2.5
cum and ~2.5-50 um)
c- air temperature above canopy

																									

c        write(5,'(22(f13.5,2x))') day_y,Soutplot1,Loutplot1,
c     .  rn4compplot1,
c     .  rnnrliteplot1,Ts4compmaizeplot1,TsIRTplot1,ea_can,
c     .  tpplot1_5,ppplot1_10,ppplot1_20,
c     .  ppplot1_30,ppplot1_40,ppplot1_60,ppplot1_100

        write(5,'(30(f13.5,2x))') day_y,Soutplot1,Loutplot1,
     .  rn4compplot1,rnnrliteplot1,NRlite680,NRlite682,NRlite698,
     .  TsIRTplot1, Ts4compmaizeplot1,	
     .  H_eddy,LE_eddy,Fco2,G1_insitu,G1_remoteavg,Tair_can,ea_can,
     .  Tset40plot1I,Tset42plot1R,Tset43plot1I,
     .  Tset44plot1I,Tset45plot1I,Tset39plot1R,Tset46plot1R,Tset38plot1R

c        write(6,'(22(f13.5,2x))') day_y,tpplot2_5,ppplot2_10,ppplot2_20,
c     .  ppplot2_30,ppplot2_40,ppplot2_60,ppplot2_100

        write(6,'(30(f13.5,2x))') day_y,Soutplot2,Loutplot2,
     .  rn4compplot2,rnnrliteplot2,NRlite680,NRlite682,NRlite698,
     .  TsIRTplot2,Ts4compmaizeplot2,	
     .  H_eddy,LE_eddy,Fco2,G2_insitu,G2_remoteavg,Tair_c2,ea_c2,
     .  Tset06plot2R,Tset08plot2I,Tset07plot2I,
     .  Tset09plot2R,Tset20plot2R,Tset23plot2I,Tset24plot2I,Tset62plot2R



      goto 1000
c        if(RH.ne.-99.99.and.Tair.ne.-99.99.and.cmup306.ne.-99.99.and.
c     .  cgup306.ne.-99.99.and.cmdn306.ne.-99.99.and.cgdn306.ne.-99.99.
c     .  and.prt306.ne.-99.99) then

          ea = es(Tair)*(RH/100.0)

C---------Total sky emissivity according to Brutsaert

	    esky0= 1.24*(ea/(Tair+273.15))**0.143
          longinfromta = esky0*sigma*(Tair+273.15)**4.0

          epsa = longin/(sigma*(Tair+273.15)**4.0)


          if((((epsa/esky0)-1.0)/0.2).gt.0.0) then
            cloudcov = sqrt(((epsa/esky0)-1.0)/0.2)
          else
            cloudcov = 0.0
          endif





c       else
        Sin = -99.99
        Sout = -99.99
        longin = -99.99
        longout = -99.99
        rnet4comp = -99.99
        rnet686 = -99.99
        longinfromta = -99.99
        esky0 = -99.99
        epsa = -99.99
        cloudcov = -99.99
c        endif
c        if(daytim.ge.227.and.daytim.le.274) then
c        write(5,'(22(f13.5,2x))') day_y,dectim,Sin,Sout,longin,longout,
c     .  rnet4comp,rnet686,longinfromta,esky0,epsa,cloudcov,Tair,RH,
c     .  Tsmaize
c        endif

1000    continue
        enddo

        STOP
        END


	  FUNCTION ES(T)

*     Deze functie berekent de verzadigingsdampspanning

	  REAL A,B,E0,ES,T
	    A = 7.5
	    B = 237.3
	    E0 = 6.107
	    ES = E0 * 10**(A*T/(B+T))
	  RETURN
      END   



	SUBROUTINE SINUSB(DAYNR1,TIMENEW,AVGGLOB,ATMTR,FRDIF,SINB,DEC)
C*********************************************
	IMPLICIT REAL(A-Z)

        integer daynr1



C-------Constants

	PI = 3.1415927
	LAT = 51.5
        SOLAR = 1370


C-------Shortwave radiation

C-------Auxiliary variables

	DAG = DAYNR1
	HOUR = TIMENEW

	RAD = PI/180

	DEC = -23.45*COS(PI*(DAG+10.0)/183)
	SINLD = SIN(RAD*LAT)*SIN(RAD*DEC)
	COSLD = COS(RAD*LAT)*COS(RAD*DEC)

C-------Sine of the solar height (inclination)

	SINB = SINLD + COSLD*COS(2*PI*(HOUR+11.0)/24.0)

        IF(SINB.LT.0) THEN
          SINB = 0.01
        ENDIF


C------CRAD is instanteneous global radiation in W m-2
C------ATMTR is instanteneous atmospheric transmissivity

        CRAD = AVGGLOB
	ATMTR = CRAD/(SOLAR*SINB)

C------Fraction diffuse

        IF(ATMTR.LE.0.22) THEN
	  FRDIF = 1.0
        ELSEIF(ATMTR.GT.0.22.AND.ATMTR.LE.0.35) THEN
          FRDIF = 1.0 - 6.4*(ATMTR-0.22)**2
        ELSE
          FRDIF = 1.47 - 1.66*ATMTR
        ENDIF


	RETURN
	END
   


