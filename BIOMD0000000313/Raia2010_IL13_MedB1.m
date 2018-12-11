function xdot = Raia2010_IL13_MedB1(time, x_values)
% function Raia2010_IL13_MedB1 takes
%
% either	1) no arguments
%       	    and returns a vector of the initial values
%
% or    	2) time - the elapsed time since the beginning of the reactions
%       	   x_values    - vector of the current values of the variables
%       	    and returns a vector of the rate of change of value of each of the variables
%
% Raia2010_IL13_MedB1 can be used with MATLABs odeN functions as 
%
%	[t,x] = ode23(@Raia2010_IL13_MedB1, [0, t_end], Raia2010_IL13_MedB1)
%
%			where  t_end is the end time of the simulation
%
%The variables in this model are related to the output vectors with the following indices
%	Index	Variable name
%	  1  	  Rec
%	  2  	  Rec_i
%	  3  	  IL13_Rec
%	  4  	  p_IL13_Rec
%	  5  	  p_IL13_Rec_i
%	  6  	  JAK2
%	  7  	  pJAK2
%	  8  	  SHP1
%	  9  	  STAT5
%	  10  	  pSTAT5
%	  11  	  SOCS3mRNA
%	  12  	  DecoyR
%	  13  	  IL13_DecoyR
%	  14  	  SOCS3
%	  15  	  CD274mRNA
%	  16  	  IL13
%
%--------------------------------------------------------
% output vector

xdot = zeros(16, 1);

%--------------------------------------------------------
% compartment values

cell = 1;

%--------------------------------------------------------
% parameter values

IL13stimulation = 1;
Kon_IL13Rec = 0.00341992;
Rec_phosphorylation = 999.631;
pRec_intern = 0.15254;
pRec_degradation = 0.172928;
Rec_intern = 0.103346;
Rec_recycle = 0.00135598;
JAK2_phosphorylation = 0.157057;
pJAK2_dephosphorylation = 0.000621906;
STAT5_phosphorylation = 0.0382596;
pSTAT5_dephosphorylation = 0.000343392;
SOCS3mRNA_production = 0.00215826;
DecoyR_binding = 0.000124391;
JAK2_p_inhibition = 0.0168268;
SOCS3_translation = 11.9086;
SOCS3_accumulation = 3.70803;
SOCS3_degradation = 0.0429186;
CD274mRNA_production = 8.21752e-05;

%--------------------------------------------------------
% initial values of variables - these may be overridden by assignment rules
% NOTE: any use of initialAssignments has been considered in calculating the initial values

if (nargin == 0)

	% initial time
	time = 0;

	% initial values
	Rec = 1.3;
	Rec_i = 113.194;
	IL13_Rec = 0;
	p_IL13_Rec = 0;
	p_IL13_Rec_i = 0;
	JAK2 = 2.8;
	pJAK2 = 0;
	SHP1 = 91;
	STAT5 = 165;
	pSTAT5 = 0;
	SOCS3mRNA = 0;
	DecoyR = 0.34;
	IL13_DecoyR = 0;
	SOCS3 = 0;
	CD274mRNA = 0;
	IL13 = 2.265;

else
	% floating variable values
	Rec = x_values(1);
	Rec_i = x_values(2);
	IL13_Rec = x_values(3);
	p_IL13_Rec = x_values(4);
	p_IL13_Rec_i = x_values(5);
	JAK2 = x_values(6);
	pJAK2 = x_values(7);
	SHP1 = x_values(8);
	STAT5 = x_values(9);
	pSTAT5 = x_values(10);
	SOCS3mRNA = x_values(11);
	DecoyR = x_values(12);
	IL13_DecoyR = x_values(13);
	SOCS3 = x_values(14);
	CD274mRNA = x_values(15);
	IL13 = x_values(16);

end;

%--------------------------------------------------------
% assignment rules
IL13 = 2.265*IL13stimulation;

%--------------------------------------------------------
% algebraic rules

%--------------------------------------------------------
% calculate concentration values

if (nargin == 0)

	% initial values
	xdot(1) = 1.3;
	xdot(2) = 113.194;
	xdot(3) = 0;
	xdot(4) = 0;
	xdot(5) = 0;
	xdot(6) = 2.8;
	xdot(7) = 0;
	xdot(8) = 91;
	xdot(9) = 165;
	xdot(10) = 0;
	xdot(11) = 0;
	xdot(12) = 0.34;
	xdot(13) = 0;
	xdot(14) = 0;
	xdot(15) = 0;
	xdot(16) = IL13;

else

	% rate equations
	xdot(1) = ( - (Kon_IL13Rec*IL13*Rec*cell) - (Rec_intern*Rec*cell) + (Rec_recycle*Rec_i*cell))/cell;
	xdot(2) = ( + (Rec_intern*Rec*cell) - (Rec_recycle*Rec_i*cell))/cell;
	xdot(3) = ( + (Kon_IL13Rec*IL13*Rec*cell) - (Rec_phosphorylation*IL13_Rec*pJAK2*cell))/cell;
	xdot(4) = ( + (Rec_phosphorylation*IL13_Rec*pJAK2*cell) - (pRec_intern*p_IL13_Rec*cell))/cell;
	xdot(5) = ( + (pRec_intern*p_IL13_Rec*cell) - (pRec_degradation*p_IL13_Rec_i*cell))/cell;
	xdot(6) = ( - (JAK2_phosphorylation*IL13_Rec*JAK2/(1+JAK2_p_inhibition*SOCS3)*cell) - (JAK2_phosphorylation*p_IL13_Rec*JAK2/(1+JAK2_p_inhibition*SOCS3)*cell) + (pJAK2_dephosphorylation*pJAK2*SHP1*cell))/cell;
	xdot(7) = ( + (JAK2_phosphorylation*IL13_Rec*JAK2/(1+JAK2_p_inhibition*SOCS3)*cell) + (JAK2_phosphorylation*p_IL13_Rec*JAK2/(1+JAK2_p_inhibition*SOCS3)*cell) - (pJAK2_dephosphorylation*pJAK2*SHP1*cell))/cell;
	xdot(8) = 0;
	xdot(9) = ( - (STAT5_phosphorylation*STAT5*pJAK2*cell) + (pSTAT5_dephosphorylation*pSTAT5*SHP1*cell))/cell;
	xdot(10) = ( + (STAT5_phosphorylation*STAT5*pJAK2*cell) - (pSTAT5_dephosphorylation*pSTAT5*SHP1*cell))/cell;
	xdot(11) = ( + (pSTAT5*SOCS3mRNA_production*cell))/cell;
	xdot(12) = ( - (DecoyR_binding*IL13*DecoyR*cell))/cell;
	xdot(13) = ( + (DecoyR_binding*IL13*DecoyR*cell))/cell;
	xdot(14) = ( + (SOCS3mRNA*SOCS3_translation/(SOCS3_accumulation+SOCS3mRNA)*cell) - (SOCS3_degradation*SOCS3*cell))/cell;
	xdot(15) = ( + (pSTAT5*CD274mRNA_production*cell))/cell;
	xdot(16) = 0;

end;
