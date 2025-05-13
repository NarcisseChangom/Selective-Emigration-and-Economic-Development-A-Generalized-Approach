

%-------------------------------------------------------------------------*
%      SELECTIVE MIGRATION, HUMAN CAPITAL AND DEVELOPMENT:                *
%                   A GENERALIZED APPROACH                                *
%                   Authors:
%                       - Narcisse Cha'ngom (LISER)     
%                       - Christoph Deuster (European Union)
%                       - Frédéric Docquier (LISER) 
%                       - Joel Machado (LISER) 
%                   Last update: January 03, 2025                         *
%-------------------------------------------------------------------------*
clc 
clear all   

global o d U dums MijH_data MijL_data 
global Nh_data Nl_data N_data Niih_data Niil_data Nii_data  Lh_data Ll_data L_data 
 
%%%
% Load raw data 
%%%
migdta = readtable('data\Movers_IER.csv','Headerlines',0,'ReadVariableNames',true);
migdta = sortrows(migdta,{'isoo','isod'}); 

locationdta = readtable('data\Location_charac_IER.csv','Headerlines',0,'ReadVariableNames',true);
locationdta = sortrows(locationdta,{'isoo'}); 

%%% 0- Identifiers 
o  = height(unique(migdta.isoo));
d  = height(unique(migdta.isod));
U= height(migdta); 
names={'Afghanistan','Angola',	'Albania','United Arab Emirates','Argentina',...
    'Armenia', 'Australia','Austria',	'Azerbaijan', 'Burundi', 'Belgium',	'Benin',...
    'Burkina Faso','Bangladesh', 'Bulgaria',	'Bahrain',	'Bahamas, The',	...
    'Bosnia and Herzegovina',	'Belarus',	'Belize',	'Bolivia',	'Brazil',...
    'Barbados',	'Brunei',	'Bhutan',	'Botswana',	'Central African Republic',	'Canada',...
    'Switzerland',	'Chile',	'China',	'Cote d`Ivoire',	'Cameroon',	'Congo, Rep. of the',...
    'Colombia',	'Comoros',	'Cape Verde',	'Costa Rica',	'Cuba',	'Cyprus',	'Czech Republic',...
    'Germany',	'Djibouti',	'Denmark',	'Dominican Republic',	'Algeria',	'Ecuador',	'Egypt',...
    'Eritrea',	'Spain',	'Estonia',	'Ethiopia',	'Finland',	'Fiji',	'France',	'Micronesia, Federated States of',...
    'Gabon',	'United Kingdom',	'Georgia',	'Ghana',	'Guinea',	'Gambia, The',	'Guinea-Bissau',...
    'Equatorial Guinea',	'Greece',	'Grenada',	'Guatemala',	'Guyana',	'Honduras',	'Croatia',	'Haiti',...
    'Hungary',	'Indonesia',	'India',	'Ireland',	'Iran',	'Iraq',	'Iceland',	'Israel',	'Italy',...
    'Jamaica',	'Jordan',	'Japan',	'Kazakhstan',	'Kenya',	'Kyrgyzstan',	'Cambodia',	'Kuwait',...
    'Laos',	'Lebanon',	'Liberia',	'Libya',	'Saint Lucia',	'Sri Lanka',	'Lesotho',	'Lithuania',...
    'Luxembourg',	'Latvia',	'Morocco',	'Moldova',	'Madagascar',	'Maldives',	'Mexico',	'Macedonia',...
    'Mali',	'Malta',	'Burma (Myanmar)',	'Mongolia',	'Mozambique',	'Mauritania',	'Mauritius',	'Malawi',...
    'Malaysia',	'Namibia',	'Niger',	'Nigeria',	'Nicaragua',	'Netherlands',	'Norway',	'Nepal',...
    'New Zealand',	'Oman',	'Pakistan',	'Panama',	'Peru',	'Philippines',	'Papua New Guinea',	'Poland',...
    'Portugal',	'Paraguay',	'Qatar',	'Romania',	'Rwanda',	'Saudi Arabia',	'Sudan',	'Senegal',...
    'Singapore',	'Solomon Islands',	'Sierra Leone',	'El Salvador',	'Somalia',	'Serbia and Montenegro',...
    'Sao Tome and Principe',	'Suriname',	'Slovakia',	'Slovenia',	'Sweden',	'Swaziland',	'Syria',...
    'Chad',	'Togo',	'Thailand',	'Tajikistan',	'Turkmenistan',	'Tonga',	'Trinidad and Tobago',...
    'Tunisia',	'Turkey',	'Tanzania',	'Uganda',	'Ukraine',	'Uruguay',	'United States',...
    'Uzbekistan',	'Saint Vincent and the Grenadines',	'Venezuela',	'Vietnam',	'Vanuatu',...
    'Samoa',	'Yemen',	'South Africa',	'Congo, Dem. Rep. of the',	'Zambia',	'Zimbabwe'}';

isoo = unique(migdta.isoo);

%% 1 - Dyadic characteristics 
dums  = reshape(migdta.dums,[o,d]); 
MijH_data  = reshape(migdta.Nijh,[o,d]);
MijL_data =  reshape(migdta.Nijl,[o,d]); 

%% 2 - Location characteristics 
Nh_data = sum(MijH_data,1);
Nl_data = sum(MijL_data,1);
N_data = Nh_data+Nl_data;
Hn_data = Nh_data./N_data;
Niih_data = sum(dums.*MijH_data,1);
Niil_data = sum(dums.*MijL_data,1);
Nii_data=Niih_data+Niil_data;
Lh_data = sum(MijH_data,2)';
Ll_data = sum(MijL_data,2)';
L_data = Lh_data+Ll_data;
Hr_data = Lh_data./L_data;
wH_data = locationdta.wH_data';
wL_data = locationdta.wL_data';
w_data = locationdta.w_data'; 
RemiT_data = locationdta.RemiT_data';
GovExp_data = locationdta.GovExp_data';
EduExp_data = locationdta.EduExp_data';
GovCons_data = locationdta.GovCons_data';
GNIc_data = locationdta.GDPc_data';
GDPc = locationdta.GDPc_data';
c_high = locationdta.c_high'; % yearly public cost of training one student in tertiary education
c_low  = locationdta.c_low'; 
CRi = c_high./c_low;


