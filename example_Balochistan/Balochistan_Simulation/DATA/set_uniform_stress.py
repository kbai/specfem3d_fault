#!/usr/python
import os
from math import *
Mega=1.0e6
sigma_1=float(raw_input("sigma_SOUTH_1?(the maximum compression(negative)(Mpa))"))
sigma_1=sigma_1*Mega
sigma_2=float(raw_input("sigma_SOUTH_2?(the other sigma in horizontal direction)(Mpa)"))
sigma_2=sigma_2*Mega
sigma_3=float(raw_input("sigma_SOUTH_3?(the sigma value in vertical direction)(Mpa)"))
sigma_3=sigma_3*Mega


azimuth=float(raw_input("the azimuth of sigma_1?(0-360)"))
az=(azimuth-45)/180*pi

sigma_yy=sigma_1*cos(az)*cos(az)+sigma_2*sin(az)*sin(az)
sigma_xx=sigma_2*cos(az)*cos(az)+sigma_1*sin(az)*sin(az)
sigma_xy=(sigma_1-sigma_2)*sin(az)*cos(az)
sigma_zz=sigma_3

settings="&UNIFORM_SOUTH Sigma_SOUTH=%f,%f,%f,%f,0.0,0.0,GradientZ = 0.0/ !stress field at SOUTH favors dip-slip motion"%(sigma_xx,sigma_yy,sigma_zz,sigma_xy)
s = "awk '{if($0~/UNIFORM_SOUTH/) printf(\"%s\\n\");else print $0}'  Par_file_faults > Par_file_faults_new" %(settings)
os.system(s)

Mega=1.0e6
sigma_1=float(raw_input("sigma_NORTH_1?(the maximum compression(negative)(Mpa))"))
sigma_1=sigma_1*Mega
sigma_2=float(raw_input("sigma_NORTH_2?(the other sigma in horizontal direction)(Mpa)"))
sigma_2=sigma_2*Mega
sigma_3=float(raw_input("sigma_NORTH_3?(the sigma value in vertical direction)(Mpa)"))
sigma_3=sigma_3*Mega


azimuth=float(raw_input("the azimuth of sigma_1?(0-360)"))
az=(azimuth-45)/180*pi

sigma_yy=sigma_1*cos(az)*cos(az)+sigma_2*sin(az)*sin(az)
sigma_xx=sigma_2*cos(az)*cos(az)+sigma_1*sin(az)*sin(az)
sigma_xy=(sigma_1-sigma_2)*sin(az)*cos(az)
sigma_zz=sigma_3

settings="&UNIFORM_NORTH Sigma_NORTH=%f,%f,%f,%f,0.0,0.0,GradientZ = 0.0/ !stress field at NORTH favors strike-slip motion"%(sigma_xx,sigma_yy,sigma_zz,sigma_xy)
s = "awk '{if($0~/UNIFORM_NORTH/) printf(\"%s\\n\");else print $0}'  Par_file_faults_new > Par_file_faults_new_2" %(settings)
os.system(s)
os.system("mv Par_file_faults_new_2 Par_file_faults")
