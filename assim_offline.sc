#!/bin/bash

if [ $# -gt 0 ] || [ $# -lt 0 ]; then
    echo "Usage:"
    echo "   ./assim_offline.sc"
    echo ""  
    echo "   0 argument needed."  
  exit
fi

# Progress bar
# 1. Create ProgressBar function
# 1.1 Input is currentState($1) and totalState($2)
function ProgressBar {
# Process data
    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
# Build progressbar string lengths
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

# 1.2 Build progressbar strings and print the ProgressBar line
# 1.2.1 Output example:                           
# 1.2.1.1 Progress : [########################################] 100%
printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%%"

}

#-------------
# Parameters |
#-------------
first_year=1958
last_year=2000 
exp_name="${first_year}-${last_year}_all_vars_1dot5std_SH_500km-grid_sx200_xy100"
moddata_co_file=moddata_co_1dot5std
make_posterior="True" # To reconstruct the posterior
outfolder_rec='/nas07/dalaiden/cyfast/paleoPF_ant/DA_exps_outputs'
frequence_sampling="3" # Frequency of the sampling. Will determine how big will the ensemble be. In season, multiple of 4.

echo "continue the work to take into the last year"
exit 1

# Full location of the folder containing the model input files (without the last "/"), and realm of the variable to be assimilated (ocean or atmos) 
here=$(pwd)
address_model_ensemble=("${here}/input/var_accu/model/files/ensemble"
                        "${here}/input/var_d18Op/model/files/ensemble")
declare -a directory_input=(""${address_model_ensemble[0]}" : atmos"
                            ""${address_model_ensemble[1]}" : atmos2") # two variables assimilated

###################################################################
# DON'T EDIT NEXT LINES
###################################################################

# Length of input files and frequency of the assimilation
duration_data="1026" # length of data to be assimilated (in seasons). Can be different than duration_model. If duration_data > duration_model, increase_ensemble_size has to be equal to "1". 
frequence_assim="1" # frequency of assim: 1 (monthly) or 12 (annual). Can be more (multiple of 12).

# Increase ensemble size options
increase_ensemble_size="1" # "1" to increase the ensemble size by selecting particles that belong to other years than the actual year of experiment, "0" to limit the ensemble size to the number of model simulations actually available.

# Area of oceanic cells (only when assimilating oceanic boxes)
fcostEPFfile=""
OceanAreaPath=""
OceanAreaVar="areacello"

module purge
module load CDO
first_file=$(find "${address_model_ensemble[0]}" -type f -name "*.nc" | head -n 1)
duration_model=`cdo ntime $first_file`
nb_simus=$(ls -l "${address_model_ensemble[0]}" | grep "^-" | wc -l)

source src/modules.load
echo "-------------------------"

echo "Experiment name: ${exp_name}"
echo "Period: ${first_year}-2025 (annual assimilation)"
echo "Namelist: ${moddata_co_file}"
nb_particles=$( expr $nb_simus '*' $duration_model '/' "$frequence_sampling" )
echo "Number of particles: $nb_particles"

offset_data=`expr ${first_year} - 1000`
offset_data=$( expr 1 '*' "$offset_data" )

rm -rf rundir/$exp_name

mkdir -p rundir/$exp_name

# Copy the information related to the prior
cp -r info_prior rundir/$exp_name/.
echo "$first_year" > rundir/$exp_name/first_year_rec

cd rundir/$exp_name

# get length of an array
arraylength=${#directory_input[@]}

# use for loop to read all values and indexes
echo 'Variables taken into account'
for (( i=1; i<${arraylength}+1; i++ ));
do
  nbparticles=0
  echo "  ${directory_input[$i-1]}"
  indir=$(echo ${directory_input[$i-1]} | cut -d : -f 1 | xargs)
  intype=$(echo ${directory_input[$i-1]} | cut -d : -f 2 | xargs)
  mkdir -p wkdir/${intype}
  ln -s ${indir}/*.nc wkdir/${intype}
  for file in wkdir/${intype}/*.nc; do let "nbparticles++"; mv "$file" "${file/.nc/_${intype}.nc}"; done
  [[ ${intype} == "atmos" ]] && liste_files_atmos=( `ls -v wkdir/${intype}/*.nc` )
  [[ ${intype} == "atmos2" ]] && liste_files_atmos2=( `ls -v wkdir/${intype}/*.nc` ) 
  [[ ${intype} == "atmos3" ]] && liste_files_atmos3=( `ls -v wkdir/${intype}/*.nc` ) 
  [[ ${intype} == "atmos4" ]] && liste_files_atmos4=( `ls -v wkdir/${intype}/*.nc` )
  [[ ${intype} == "atmos5" ]] && liste_files_atmos5=( `ls -v wkdir/${intype}/*.nc` )
  [[ ${intype} == "atmos6" ]] && liste_files_atmos6=( `ls -v wkdir/${intype}/*.nc` ) 
  [[ ${intype} == "ocean" ]] && liste_files_ocean=( `ls -v wkdir/${intype}/*.nc` )
  [[ ${intype} == "evolu" ]] && liste_files_evolu=( `ls -v wkdir/${intype}/*.nc` )
done
echo "-------------------------"

idx_number_assim=$((offset_data/frequence_assim))

progress_idx=0
end_progress_bar=$((duration_data-offset_data))
for timestep in $(seq $((1+$offset_data)) $frequence_assim $duration_data)
do
        # Progress bar
        ProgressBar $progress_idx $end_progress_bar
        
        idx_num_part=1
 
        month_start_data=$(($timestep))
        month_end_data=$(($timestep+$frequence_assim-1))

        if (( "$timestep" <= "$duration_model" )); then
            
            month_start_model=$month_start_data
            month_end_model=$month_end_data
            
            for (( idx=1; idx<${nbparticles}+1; idx++ )); # boucle sur les particules
            do

                file_name_atmos=${liste_files_atmos[$idx-1]}
                file_name_atmos2=${liste_files_atmos2[$idx-1]}
                file_name_atmos3=${liste_files_atmos3[$idx-1]}
                file_name_atmos4=${liste_files_atmos4[$idx-1]}
                file_name_atmos5=${liste_files_atmos5[$idx-1]}
                file_name_atmos6=${liste_files_atmos6[$idx-1]}
                file_name_ocean=${liste_files_ocean[$idx-1]}
                file_name_evolu=${liste_files_evolu[$idx-1]}
                        
                sed "s,{file_newdata_atmos},${file_name_atmos},g;\
                     s,{file_newdata_atmos2},${file_name_atmos2},g;\
                     s,{file_newdata_atmos3},${file_name_atmos3},g;\
                     s,{file_newdata_atmos4},${file_name_atmos4},g;\
                     s,{file_newdata_atmos5},${file_name_atmos5},g;\
                     s,{file_newdata_atmos6},${file_name_atmos6},g;\
                     s,{file_newdata_ocean},${file_name_ocean},g;\
                     s,{file_newdata_evolu},${file_name_evolu},g;\
                     s,{NTIMERUNY},1,g;\
                     s,{NTIMERUND},1,g;\
                     s,{FCOSTEPFOUT},${fcostEPFfile},g;\
                     s,{CLIOareaPath},${OceanAreaPath},g;\
                     s,{CLIOareaVar},${OceanAreaVar},g;\
                     s,{frequence_mod},${frequence_assim},g;\
                     s,{nb_variables},${arraylength},g;\
                     s,{debut_mois_data},${month_start_data},g;\
                     s,{fin_mois_data},${month_end_data},g;\
                     s,{debut_mois_model},${month_start_model},g;\
                     s,{fin_mois_model},${month_end_model},g;" < ../../$moddata_co_file > moddata_co.namelist

                mkdir -p output_log
                ../../src/moddata_co FAP 1 off >> output_log/moddata_co_$idx.log
                
                save_name=${idx_num_part}
                mkdir -p fcost_folder_tmp/$save_name
                mv fcostopt.dat fcost_folder_tmp/$save_name/
                mv mois_start.dat fcost_folder_tmp/$save_name/
                mv mois_end.dat fcost_folder_tmp/$save_name/
                
                idx_num_part=$((idx_num_part+1))               
            done
        fi

    	month_start_model=1
    	month_end_model=$(($month_start_model+$frequence_assim-1))
  
        for idate in $(seq 1 1 $((((duration_model*increase_ensemble_size)/frequence_sampling))))
        do
                
            if (( month_start_model == month_start_data )) && (( "$timestep" <= "$duration_model" )); then # on ne va pas repiocher ds l'année en cours, sauf qd au est au delà de la période couverte par le modèle
                var_tmp='aaa'
            else
                for (( idx=1; idx<${nbparticles}+1; idx++ )); # boucle sur les particules
                do                    

                    file_name_atmos=${liste_files_atmos[$idx-1]}
                    file_name_atmos2=${liste_files_atmos2[$idx-1]}
                    file_name_atmos3=${liste_files_atmos3[$idx-1]}
                    file_name_atmos4=${liste_files_atmos4[$idx-1]}
                    file_name_atmos5=${liste_files_atmos5[$idx-1]}
                    file_name_atmos6=${liste_files_atmos6[$idx-1]}
                    file_name_ocean=${liste_files_ocean[$idx-1]}
                    file_name_evolu=${liste_files_evolu[$idx-1]}

                    sed "s,{file_newdata_atmos},${file_name_atmos},g;\
                            s,{file_newdata_atmos2},${file_name_atmos2},g;\
                            s,{file_newdata_atmos3},${file_name_atmos3},g;\
                            s,{file_newdata_atmos4},${file_name_atmos4},g;\
                            s,{file_newdata_atmos5},${file_name_atmos5},g;\
                            s,{file_newdata_atmos6},${file_name_atmos6},g;\
                            s,{file_newdata_ocean},${file_name_ocean},g;\
                            s,{file_newdata_evolu},${file_name_evolu},g;\
                            s,{NTIMERUNY},1,g;\
                            s,{NTIMERUND},1,g;\
                            s,{FCOSTEPFOUT},${fcostEPFfile},g;\
                            s,{CLIOareaPath},${OceanAreaPath},g;\
                            s,{CLIOareaVar},${OceanAreaVar},g;\
                            s,{frequence_mod},${frequence_assim},g;\
                            s,{nb_variables},${arraylength},g;\
                            s,{debut_mois_data},${month_start_data},g;\
                            s,{fin_mois_data},${month_end_data},g;\
                            s,{debut_mois_model},${month_start_model},g;\
                            s,{fin_mois_model},${month_end_model},g;" < ../../$moddata_co_file > moddata_co.namelist
        
                    mkdir -p output_log
                    ../../src/moddata_co FAP 1 off >> output_log/moddata_co_$idx.log
                                                    
                    save_name=${idx_num_part}
                    mkdir -p fcost_folder_tmp/$save_name
                    mv fcostopt.dat fcost_folder_tmp/$save_name/
                    mv mois_start.dat fcost_folder_tmp/$save_name/
                    mv mois_end.dat fcost_folder_tmp/$save_name/	
                    
                    idx_num_part=$((idx_num_part+1))               
                done     
                
            fi
            month_start_model=$((month_start_model+frequence_sampling))
            month_end_model=$((month_end_model+frequence_sampling))  
                
        done #idate
        
        nbr_fcost=$(find fcost_folder_tmp/* -maxdepth 0 -type d | wc -l)
        mkdir -p output_fcost
        ../../src/PartFilter `pwd`/fcost_folder_tmp/ . $nbr_fcost FAP -1 > output_fcost/fcost_month_${timestep}
        
        # on change le nom en fonctionque l'assim soit annual ou mensuelle
        if (( $frequence_assim == 12 )); then
            mv output_fcost/fcost_month_${timestep} output_fcost/fcost_year_${idx_number_assim}
        fi

        idx_number_assim=$((idx_number_assim+1))

        # Progress bar
        progress_idx=$((progress_idx+1))
done

ProgressBar $end_progress_bar $end_progress_bar
echo

if [ "$make_posterior" == "True" ]; then

    echo '-----------------------'
    echo "Reconstruct the fields|"
    echo '-----------------------'

    cd ${here}/post
    export PYTHONWARNINGS="ignore"
    module purge
    module load releases/2022b
    module load ELIC_Python/1-foss-2022b
    python -W ignore make_prior.py $exp_name

    python -W ignore make_posterior.py $exp_name $outfolder_rec

    cd ${here}/rundir/$exp_name

fi

/bin/rm -rf fcost_folder_tmp/*
/bin/rm -rf wkdir
