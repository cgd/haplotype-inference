#!/bin/bash
true_script_location=`readlink -fn $0`

pushd `dirname $true_script_location`

for  i in `"ls" lib`;
do
  CLASSPATH=$CLASSPATH:lib/$i;
done

java -cp $CLASSPATH:dist/haplotype-inference-1.0.jar org.jax.haplotype.io.SnpStreamUtil

popd
