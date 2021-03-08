molecule="TMETN"


echo "## Step 1 : SMILES to PDB"
## Using openbabel
#export PATH=${PATH}:"/g/g92/pham20/codes/install_openbabel/"
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/g/g92/pham20/codes/install_openbabel/lib64"
obabel --gen3D -ismiles  ${molecule}.smiles -o pdb &>   ${molecule}.pdb 


echo "## Step 2 : Generate the *.key file"
echo "parameters   $HOME/codes/tinker/params/mm3" > ${molecule}.key
#echo "a-axis                    40.00" >> ${molecule}.key 
#echo "b-axis                    40.00" >> ${molecule}.key 
#echo "c-axis                    40.00" >> ${molecule}.key 



echo "## Step 3 : pdb to TINKER XYZ; ready to run"
$HOME/codes/tinker/source/pdbxyz.x     ${molecule}.pdb    -k ${molecule}.key
