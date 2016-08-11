function ps = inx(vr,Z);

[nv,dm]=size(Z);

i=1;
ps=0;

while i<=nv
   
  tmp=Z{i};  
  tst1=strfind(tmp,'(');
  tst2=isempty(tst1);
  if tst2==0
      [vrnm,per]=strtok(tmp,'(');
  else
      vrnm=deblank(tmp);
  end   
  cmp=strcmp(vr,vrnm);
  
  if cmp==0
      i=i+1;
  else
      ps=i;
      i=nv+1;
  end
  
end  
  