pro compute_Fig_S7_S8

  ;.r ~/Lib/astron/pro/legend.pro


  ; NOTE; MDOT IS DEFINED NEGATIVE
  
  grav = 6.67259d-8 ; [cm3 g-1 s-2]
  M_SUN    =1.989d+33 ; g
  AU       =1.496d13 ; cm
  ysec     =31557600d0 ; sec
  R_sun    =6.96d10 ; cm
  L_SUN    =3.826d+33 ; erg s-1
  c        =2.9979245800d10; cm/s


  ; initialize some output parameters
  change_rad_min = 1d5
  change_rad_max = 0.

  ; take parameters as Chen et al. 2019
  R_star = 267.*R_sun ; in cm

 
  ; beta velocity profile
  beta_all = [0.5, 5.] ; C-rich, O-rich

  ; ADAPT INPUT
  ; mass of primary
  M1_Msun = 1.5 ;Msun 1.5 or 3.   [should be larger than 1.2 Msun]
  M1 = M1_Msun * M_sun
  M_final = 0.6 ; see Cummings 2018   0.6 or 0.7; depending on initial mass
  input_scale=1.25 ; or 0.0179+0.733*M1_Msun  ;  max widening is 0.0179+0.733*M1_Msun  , max shrinking is a_ini/R_star ; 1.5 or 5.

  ;set_plot,'ps'
  ;device,/color,file='Fig_S7.eps',/encapsulated,yoff=0,/portrait, xoff=0,xsize=15,ysize=15, set_font='Courier',/tt_font
  

  ; mass of secondary
  m2_msun = [1.2, 0.6, 0.3, 0.01] ; Msun
  index = ['(A)','(B)','(C)','(D)']
  ;m2Msun = [1.2, 0.3, 0.001, 3e-6]; Msun
  m2 = m2_msun * M_sun


  ; initial orbital separation
  a_ini_all = [2., 4., 10., 25.]
  scale_plot=findgen(n_elements(a_ini_all))
  for is=0,n_elements(scale_plot)-1 do scale_plot(is)=input_scale


  ; initial q
  q = M1/m2
  q_string = strcompress(string(q,format='(f7.1)'),/remove_all)
  q_float = float(q_string)
  ;print, q_float, fix(q_float*10) mod 10
  for iq=0,n_elements(q_string)-1 do begin
    if ((fix(q_float(iq)*10) mod 10) eq 0) then begin
      q_string(iq) = strcompress(string(q_float(iq),format='(I3)'),/remove_all)
    endif else begin
      q_string(iq) = strcompress(string(q(iq),format='(f7.1)'),/remove_all)
    endelse    
  endfor
  for iq=0,n_elements(q_string)-1 do q_string(iq) = "q'="+q_string(iq)
  colors_q = ['black','goldenrod','blue','magenta']
  
  ;
   ; INITIALIZE MDOT
  ; calculate following Vassiliadis and Wood (1993) + McDonald 2016,2018
  ;assume vexp (=v_infty) and teff
  Teff = 2500 ; K
  v_infty = 15d5 ; cm/s
  logP = -2.07 + 1.94*alog10(R_star/R_sun) -0.9*alog10(M1/M_sun) ; in days
  P = 10.^(logP) ; days
  print, 'initial P is ', P, ' days'
  Mbol = -3.00 * logP + 2.85 ; see De Beck et al. 2010 from Feast et al. 1989
  logL_Lsun = (Mbol -4.75)/(-2.5)
  L_star = 10.^(logL_Lsun) * L_sun
  print, 'initial L_star is ', L_star/L_sun, ' Lsun'
;  if (P ge 300) then begin
    logMdot = -11.4 + 0.0123*P ; in Msun/yr ; for period > 300 days --> minimum Mdot around 1 - 2 e-8 Msun/yr
    Mdot_Msunyr_1 = 10.^logMdot ; in Msun/yr
    Mdot1 = Mdot_Msunyr_1 * M_sun / ysec
    Mdot_single_scattering = L_star/(c * v_infty) ; in g/s
    Mdot = -1d0*min([Mdot_single_scattering, Mdot1]) ; in g/s
  ; for 60 days< period < 300 days, follow results of Iain McDonald (2016, 2018): Mdot tussen few times 1e-9 and 1.5e-7 Msun/yr
  if (P lt 300) then Mdot = -1d0 * 2e-8*M_sun/ysec  ; see Fig 4 of McDonald 2018 (+to connect to relation of VW93)
  if (P lt 60) then begin
    print, 'period < 60 days'
    STOP
  endif
  Mdot_Msunyr = Mdot/M_sun * ysec ; in Msun/yr
  Mdot_start = Mdot ; in g/s
  print, 'start at ', L_star/L_sun, P, Mdot_Msunyr

 

  M1_before_PN_Msun = M_final ; Msun ; see Cummings (2018)
  M1_before_PN = M1_before_PN_Msun  * M_sun
  total_time = (M1_Msun - M1_before_PN_Msun)/abs(Mdot_Msunyr_1) ; in yr  -- really take the extreme so also for the too low values of mdot_1 to have enough timesteps
  


  ; take 10000 time steps
  total_steps = 100000
  delta_time = ceil(total_time / total_steps) ; in yr
  delta_time_sec = delta_time * ysec ; in sec
  total_steps = total_steps + 1
  time = delta_time * findgen(total_steps) ; in yr
  time_sec = time * ysec
  
 
  a_final=fltarr(total_steps)
  R_RL1 = fltarr(total_steps) ; Roche Lobe of primary star
  Mdot_evol = fltarr(total_steps) ; change of Mdot through the evolution
  Mdot_Msunyr_evol = fltarr(total_steps) ; change of Mdot through the evolution
  M1_evol = fltarr(total_steps) ; change of M1 through the evolution
  m2_evol = fltarr(total_steps) ; change of M2 through the evolution
  L_star_evol = fltarr(total_steps) ; change of L_star through the evolution
  P_evol = fltarr(total_steps) ; change of P through the evolution

 
  
  !p.multi = [0,2,2]
  for j = 0,n_elements(a_ini_all) -1 do begin

    a_ini_au = a_ini_all(j)
    a_ini = a_ini_au * au ; in cm

    ; initialize
    a_widen_max = 1.
    
    for i = 0,n_elements(m2)-1 do begin
      
      M2_start = m2(i)

    ; initialize
      a_fraction_widen = 1.
      a_fraction_shrink = 1 
   
      
 

      for k = 0,n_elements(beta_all)-1 do begin
        beta_w_here = beta_all(k)
 
        ; initialize
        a_here = a_ini ; in cm
        M1_this_time = M1 ; in g
        M2_now = M2_start ; in g
        end_time = max(time)
        t_end = total_steps-1
        t_merge = total_steps-1
        M1_end = M1
        m2_end = m2_start
        Mdot_end = Mdot_start ; in g/s
        
 
        t=0  
        WHILE (M1_this_time/M_sun ge M1_before_PN_Msun) DO BEGIN
          

          if (t eq 0) then begin
            M1_this_time = M1
            m2_now = m2_start
            q_now = M1_this_time / m2_now
            r_RL1(0) = 0.49 * q_now^(2./3.)/(0.6*q_now^(2./3.)+alog(1.+q_now^(1./3.))) * a_here ; in cm
            a_final(0) = a_here ; in cm
            Mdot_evol(0) = Mdot_start ; in g/s
            Mdot_Msunyr_evol(0) = Mdot_start/M_sun*ysec ; in Msun/yr
            P_evol(0) = P
            L_star_evol(0) = L_star
          endif

          if (t gt 0) then begin
            logP = -2.07 + 1.94*alog10(R_star/R_sun) -0.9*alog10(M1_this_time/M_sun) ; in days
            P = 10.^(logP) ; days
            Mbol = -3.00 * logP + 2.85 ; see De Beck et al. 2010 from Feast et al. 1989
            logL_Lsun = (Mbol -4.75)/(-2.5)
            L_star = 10.^(logL_Lsun) * L_sun
            L_star_evol(t) = L_star/L_sun ; in Lsun
            P_evol(t) = P ; days            
            if (P ge 300) then begin
              logMdot = -11.4 + 0.0123*P ; in Msun/yr  for P > 300 days
              Mdot_Msunyr_1 = 10.^logMdot ; in Msun/yr
              Mdot1 = Mdot_Msunyr_1 * M_sun / ysec
              Mdot_single_scattering = L_star/(c * v_infty) ; in g/s
              Mdot_evol(t) = -1d0*min([Mdot_single_scattering, Mdot1]) ; in g/s
            endif else begin 
              Mdot_evol(t) = -1d0 * 2e-8*M_sun/ysec  ; see note Fig. 4 of McDonald 2018
            endelse
            if (P lt 60) then begin
              print, 'period < 100 days'
              STOP
            endif
            Mdot_Msunyr_evol(t) = Mdot_evol(t)/M_sun * ysec ; in Msun/yr
            a_final(t) = a_final(t-1) ; will be adapted later on but needed now for some first calculations of some parameters such as v_orb
          endif
 


          ; orbital velocity
          v_orb = sqrt(grav*(M1_this_time+m2_now)/a_final(t)) ; in cm/s at initial a

          ; beta velocity profile
          beta_w = beta_w_here
          v_0 = 2d5 ; cm/s
          ; v_wind at a_here
          v_wind = v_0 + (v_infty - v_0)*(1.- R_star/a_here)^beta_w ; in cm/s

          ; equations from Saladino et al. 2018
          c1 = max([q_now, 0.6*q_now^1.7])
          c2 = 1.5+0.3*q_now
          eta_iso = 1./((1+q_now)^2.)
          eta = min([1./(c1+(c2*v_wind/v_orb)^3.) + eta_iso,0.6])

          k1 = 1.7 + 0.3*q_now
          k2 = 0.5 + 0.2*q_now
          alpha = 0.75 + 1./(k1+(k2*v_wind/v_orb)^5.)
          alpha_bhl = 1.
          beta_bhl = alpha_bhl/(1.+q_now)^2. * v_orb^4./(v_wind*(v_wind^2. + v_orb^2.)^(3./2.))
          beta_max = min([0.3, 1.4/q_now^2.])
          beta = min([alpha * beta_bhl, beta_max])

          ; compute new M1 and m2
          M1_this_time = M1_this_time - abs(Mdot_evol(t)) * delta_time_sec ; in g
          m2_now = m2_now + beta * abs(Mdot_evol(t)) * delta_time_sec
          M1_evol(t) = M1_this_time
          M2_evol(t) = M2_now
          q_now = M1_this_time / m2_now
     
          f = 1. - beta*q_now - eta*(1.-beta)*(1.+q(i)) - (1.-beta)*q_now/(2.*(1.+q_now))

          if (t gt 0) then begin
            adot = a_final(t-1) * (-2.) * Mdot/M1_this_time * f
            a_final(t) = a_final(t-1)+ adot*delta_time_sec
            if (a_final(t) le R_star) then begin
              t_merge = min([t,t_merge])
            endif
          endif
          r_RL1(t) = 0.49 * q_now^(2./3.)/(0.6*q_now^(2./3.)+alog(1.+q_now^(1./3.))) * a_final(t) ; in cm
          
          end_time = time(t)
          M1_end = M1_evol(t)
          M2_end = M2_evol(t)
          Mdot_end = Mdot_evol(t)
          P_end = P_evol(t)
          L_star_end = L_star_evol(t)
          t_end = t

          t = t+1
          endwhile

        ; go back one step to know all quantities g
        t_end = t-2
        end_time = time(t_end)
        
        
        fraction_a = a_final(t_end)/a_final(0)
        if (fraction_a lt 1.) then a_fraction_shrink = fraction_a
        if (fraction_a gt 1.) then a_fraction_widen = fraction_a
        if (t_merge lt total_steps-1) then begin
           a_fraction_shrink = R_star/a_final(0)
        endif 
        

        change_rad = sqrt(M2_end/M2_start)-1.
        if (change_rad gt change_rad_max) then change_rad_max = change_rad
        if (change_rad lt change_rad_min) then change_rad_min = change_rad

        
        a_final_au = a_final / au
        
        if (a_fraction_widen gt a_widen_max) then a_widen_max = a_fraction_widen

        
        
        if ((i eq 0) and (k eq 0)) then begin
          t_selec = where(time le min([end_time,time(t_merge)]))
          ymin = a_ini_au/scale_plot(j)
          ymax = min([a_ini_au*scale_plot(j),a_ini_au*2.5])
          if ((R_star/au gt ymin) or (max(R_RL1(t_selec))/au gt ymin)) then ymin=0
          plot, time, a_final_au, yra=[ymin,ymax], xtitle = 'Time on AGB [yr]', ytitle = 'Orbital separation [au]',/nodata, charthick=3,/xs,/ys,xra=[0,end_time],xcharsize=0.7,ycharsize=0.9
          oplot, time(t_selec), a_final_au(t_selec), color=fsc_color(colors_q(i)),thick=3
          ;legend, index(j), position=[end_time*0.7, ymax*0.97],box=0,charthick=3,charsize=0.7
          legend, index(j), position=[end_time*0.75, ymax*0.982],box=0,charthick=3,charsize=0.7
          if ((R_star/au gt ymin) or (max(R_RL1(t_selec))/au gt ymin)) then begin
            oplot, time(0:t_end), R_RL1(0:t_end)/au, thick=3, color=fsc_color('gray')
            hor, R_star/au, color=fsc_color('gray')
            polyfill, [0, 0,end_time, end_time],[ymin, R_star/au, R_star/au, ymin], /data, orientation = 45, color=fsc_color('gray')           
          endif
          if (j eq n_elements(a_ini_all)-1) then legend,q_string,color=fsc_color(colors_q), textcolor=fsc_color(colors_q), box=0, /left, /top, charsize=0.8,charthick=3
          if (j eq n_elements(a_ini_all)-1) then legend,['C-rich','O-rich'], line=[0,2], box=0, /right, /bottom, charsize=0.7,charthick=3
        endif else begin
          t_selec = where(time le min([end_time,time(t_merge)]))
          if (k eq 0) then  oplot, time(t_selec), a_final_au(t_selec), color=fsc_color(colors_q(i)),thick=3
          if (k eq 1) then oplot, time(t_selec), a_final_au(t_selec), color=fsc_color(colors_q(i)), line=2,thick=3
        endelse

      endfor

 
    endfor

 
  endfor
  
  t_selec = where(time le end_time)
  
  
  !p.multi=0

     ;device,/close
     ;set_plot,'x'


end