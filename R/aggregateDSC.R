aggregateDSC=function(modlist){
 mt=0;
 for(j in 1:length(modlist)){
  mt=mt+modlist[[j]]$predictions
 }
 t(apply(mt,1,function(x){x/sum(x)}))
}