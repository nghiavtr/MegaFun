#!/bin/bash
istr="/path/to/MegaFun"
ostr="$PWD"

chmod 777 MegaFun.sh
chmod 777 bin/*
eval "sed -i -e 's#"$istr"#"$ostr"#g'  MegaFun.sh"
eval "sed -i -e 's#"$istr"#"$ostr"#g'  *.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g'  buildAnno/*.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g'  buildAnno/*.sh"
echo "Done!"

