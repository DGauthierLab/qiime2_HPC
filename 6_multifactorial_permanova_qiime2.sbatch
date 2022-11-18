##Adonis script: allows multifactorial PERMANOVA

crun qiime diversity adonis \
  --i-distance-matrix ${COREMETRICS}/${DMBETA}_distance_matrix${FILTTAX}${FILTCOV}${FILTMETA}.qza \
  --m-metadata-file $METAPATH \
  --p-formula ${PERMFORM} \
  --p-permutations 999 \
  --o-visualization ${COREMETRICS}/${DMBETA}-significance_adonis_${PERMFORM}_${FILTTAX}${FILTCOV}${FILTMETA}.qzv 
