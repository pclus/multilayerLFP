# Check bipolars ----------------------------------------------
 t0=100
 tf=110
 t,bip=read_channel(100,t0,tf,"bipolar_pre");
 t,ch1=read_channel(100,t0,tf,"pre");
 t,ch2=read_channel(104,t0,tf,"pre");
 t,ch3=read_channel(108,t0,tf,"pre");

 plot(bip)
 plot!((ch1-ch3)/(2*20.0),lw=0.5)

 sum(abs.(bip-(ch1-ch3)/(2*20.0)))

 # Check CSD ---------------------------------------------------
  # Check CSD ---------------------------------------------------
 t,cs=read_channel(100,t0,tf,"csd_pre");
 t,ch5=read_channel(201,t0,tf,"pre");
 t,ch3=read_channel(199,t0,tf,"pre");
 t,ch4=read_channel(200,t0,tf,"pre");
 t,ch7=read_channel(203,t0,tf,"pre");
 t,ch8=read_channel(204,t0,tf,"pre");

 plot(cs,lw=2)
 plot!((4*ch5-ch3-ch4-ch7-ch8)/(25.61^2),lw=1,xlim=(0,100))

 sum(abs.(cs-(4*ch5-ch3-ch4-ch7-ch8)/(25.61^2)))

 # distance analysis, to validate distance between sites -------
 diff=zeros(96,6);
 for i in 1:4:384
     t,ch1=read_channel(i+0,t0,tf,fl);
     t,ch2=read_channel(i+1,t0,tf,fl);
     t,ch3=read_channel(i+2,t0,tf,fl);
     t,ch4=read_channel(i+3,t0,tf,fl);
     diff[Int(floor(i/4))+1,1]=norm(ch1-ch2)
     diff[Int(floor(i/4))+1,2]=norm(ch1-ch3)
     diff[Int(floor(i/4))+1,3]=norm(ch1-ch4)
     diff[Int(floor(i/4))+1,4]=norm(ch2-ch3)
     diff[Int(floor(i/4))+1,5]=norm(ch2-ch4)
     diff[Int(floor(i/4))+1,6]=norm(ch3-ch4)
 end

