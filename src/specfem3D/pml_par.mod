	  �>  W   k820309              12.1        Ȁ�U                                                                                                           
       pml_par.f90 PML_PAR                                                     
       CUSTOM_REAL                                                                                                     4                                                                                                                                       &                                                                                                                         &                                                                                                                         &                                                                                                                         &                                                                                                                         &                                                                                       	     	                                                   
     	                                                        	                                                                    	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                   &                   &                                                                                                        	                &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                        	                &                   &                   &                                                                                                         	                &                   &                   &                                                                                     !                   	                &                   &                   &                                                                                     "                   	                &                   &                   &                                                                                     #                   	                &                   &                   &                                                                                     $                   	                &                   &                   &                                                                                     %                   	                &                   &                   &                                                                                     &                   	                &                   &                   &                                                                                     '                   	                &                   &                   &                                                                                     (                   	                &                   &                   &                                                                                     )                   	                &                   &                   &                                                                                     *                   	                &                   &                   &                                                                                     +                   	                &                   &                   &                                                                                     ,                   	                &                   &                   &                                                                                     -                   	                &                   &                   &                                                                                     .                   	                &                   &                   &                   &                   &                                                                                     /                   	                &                   &                   &                   &                   &                                                                                     0                   	                &                   &                   &                   &                   &                                                                                     1                   	                &                   &                   &                   &                                                                                     2                   	                &                   &                   &                   &                                                                                     3                   	                &                   &                   &                   &                                                                                     4                   	                &                   &                   &                   &                                                                                     5                   	                &                   &                   &                   &                                                                                     6                   	                &                   &                   &                   &                                                                                     7                   	                &                   &                   &                   &                   &                                                                                     8                   	                &                   &                   &                   &                   &                                                                                     9                   	                &                   &                   &                   &                   &                                                                                     :                   	                &                   &                   &                   &                                                                                     ;                   	                &                   &                   &                   &                                                                                     <                   	                &                   &                   &                   &                                                                                     =                   	                &                   &                   &                   &                                                                                     >                   	                &                   &                   &                   &                                                                                     ?                   	                &                   &                   &                   &                                                                                     @                   	                &                   &                   &                   &                   &                                                                                     A                   	                &                   &                   &                   &                   &                                                                                     B                   	                &                   &                   &                   &                   &                                                                                     C                   	                &                                                                                     D                   	                &                                                                                     E                   	                &                   &                   &                   &                   &                                                                                     F                   	                &                   &                   &                   &                   &                                                                                     G                   	                &                   &                   &                   &                   &                                                                                     H                   	                &                   &                   &                   &                   &                   &                                                                                     I                   	                &                   &                   &                   &                   &                                                                                     J                   	                &                   &                   &                   &                                                                                     K                   	                &                   &                   &                                                                                     L                   	                &                   &                   &                   &                   &                   &                                                                                     M                   	                &                   &                   &                   &                   &                                                                                     N                   	                &                   &                   &                   &                   &                                                                                        O                                                         P                                                       Q                                   &                                                                                      R                                   &                                                                                        S                                                         T                                                      U                   	                &                   &                                                                                     V                   	                &                   &                                              �         fn#fn    �   L   J  CONSTANTS &     q       CUSTOM_REAL+CONSTANTS    y  @       NSPEC_CPML    �  �       CPML_TO_SPEC    E  �       SPEC_TO_CPML    �  �       CPML_REGIONS    ]  �       CPML_TYPE    �  �       IS_CPML    u  @       CPML_WIDTH_X    �  @       CPML_WIDTH_Y    �  @       CPML_WIDTH_Z    5  �       D_STORE_X    	  �       D_STORE_Y    �  �       D_STORE_Z    �  �       K_STORE_X    �  �       K_STORE_Y    Y	  �       K_STORE_Z    -
  �       ALPHA_STORE_X      �       ALPHA_STORE_Y    �  �       ALPHA_STORE_Z    �  �       DISPL_OLD    M  �       PML_DUX_DXL    	  �       PML_DUX_DYL    �  �       PML_DUX_DZL    �  �       PML_DUY_DXL    =  �       PML_DUY_DYL    �  �       PML_DUY_DZL    �  �       PML_DUZ_DXL    q  �       PML_DUZ_DYL    -  �       PML_DUZ_DZL     �  �       PML_DUX_DXL_OLD     �  �       PML_DUX_DYL_OLD     a  �       PML_DUX_DZL_OLD       �       PML_DUY_DXL_OLD     �  �       PML_DUY_DYL_OLD     �  �       PML_DUY_DZL_OLD     Q  �       PML_DUZ_DXL_OLD       �       PML_DUZ_DYL_OLD     �  �       PML_DUZ_DZL_OLD #   �  �       PML_DPOTENTIAL_DXL #   A  �       PML_DPOTENTIAL_DYL #   �  �       PML_DPOTENTIAL_DZL '   �  �       PML_DPOTENTIAL_DXL_OLD '   u  �       PML_DPOTENTIAL_DYL_OLD '   1  �       PML_DPOTENTIAL_DZL_OLD "   �  �       RMEMORY_DUX_DXL_X "   �  �       RMEMORY_DUX_DYL_X "   �   �       RMEMORY_DUX_DZL_X "   �!  �       RMEMORY_DUY_DXL_X "   �"  �       RMEMORY_DUY_DYL_X "   Y#  �       RMEMORY_DUZ_DXL_X "   -$  �       RMEMORY_DUZ_DZL_X "   %  �       RMEMORY_DUX_DXL_Y "   �%  �       RMEMORY_DUX_DYL_Y "   �&  �       RMEMORY_DUY_DXL_Y "   �'  �       RMEMORY_DUY_DYL_Y "   �(  �       RMEMORY_DUY_DZL_Y "   m)  �       RMEMORY_DUZ_DZL_Y "   A*  �       RMEMORY_DUZ_DYL_Y "   +  �       RMEMORY_DUX_DXL_Z "   �+  �       RMEMORY_DUX_DZL_Z "   �,  �       RMEMORY_DUY_DYL_Z "   �-  �       RMEMORY_DUY_DZL_Z "   e.  �       RMEMORY_DUZ_DXL_Z "   Q/  �       RMEMORY_DUZ_DYL_Z "   =0  �       RMEMORY_DUZ_DZL_Z '   )1  �       POTENTIAL_ACOUSTIC_OLD /   �1  �       POTENTIAL_DOT_DOT_ACOUSTIC_OLD '   A2  �       RMEMORY_DPOTENTIAL_DXL '   -3  �       RMEMORY_DPOTENTIAL_DYL '   4  �       RMEMORY_DPOTENTIAL_DZL &   5        RMEMORY_DISPL_ELASTIC +   	6  �       RMEMORY_POTENTIAL_ACOUSTIC #   �6  �       ACCEL_ELASTIC_CPML 0   �7  �       POTENTIAL_DOT_DOT_ACOUSTIC_CPML -   �8        RMEMORY_COUPLING_AC_EL_DISPL 1   �9  �       RMEMORY_COUPLING_EL_AC_POTENTIAL 9   u:  �       RMEMORY_COUPLING_EL_AC_POTENTIAL_DOT_DOT -   a;  @       NGLOB_INTERFACE_PML_ACOUSTIC ,   �;  @       NGLOB_INTERFACE_PML_ELASTIC .   �;  �       POINTS_INTERFACE_PML_ACOUSTIC -   m<  �       POINTS_INTERFACE_PML_ELASTIC #   �<  @       B_RECLEN_PML_FIELD '   9=  @       B_RECLEN_PML_POTENTIAL    y=  �       B_PML_FIELD     >  �       B_PML_POTENTIAL 