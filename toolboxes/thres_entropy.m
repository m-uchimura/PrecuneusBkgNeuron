function Entropy_thre=thres_entropy(Entropy,P_CHIV,P_CHIV_dotsh,I_SHUFF,p_thre)
Entropy_thre_tmp(:,:,1:2)=(P_CHIV(:,:,1:2)<=p_thre).*(P_CHIV_dotsh(:,:,1:2)<=p_thre).*Entropy(:,:,1:2);
Entropy_thre_tmp(:,:,3)=(P_CHIV(:,:,3)<=p_thre).*Entropy(:,:,3);
for i=1:size(Entropy_thre_tmp,1)
Entropy_thre(i,:,1)=Suc5(Entropy_thre_tmp(i,:,1));
Entropy_thre(i,:,2)=Suc5(Entropy_thre_tmp(i,:,2));
Entropy_thre(i,:,3)=Suc5(Entropy_thre_tmp(i,:,3));
end
ent_shuff_mean=squeeze(nanmean(I_SHUFF,2));
entropy_shuffled=ent_shuff_mean(:,:,1).*((Entropy_thre(:,:,1))>0).*(Entropy_thre(:,:,2)>0);
entropy_thre_minus_spread=squeeze(Entropy_thre(:,:,1)-entropy_shuffled);
entropy_thre_minus_spread(entropy_thre_minus_spread<0)=0;
Entropy_thre(:,:,1)=entropy_thre_minus_spread;