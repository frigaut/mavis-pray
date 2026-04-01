/*
// Actual optics configuration of mavis:
// from https://docs.google.com/spreadsheets/d/1xxNe21pZ5DOW01jhxpFV5pQQLP1y0N2zne8CT-KiV0c/edit?gid=0#gid=0
// # optics_name    altitude  nmerror checked comment
start
 1  Field_lens	      45.6    10.0
 2  DM_High           13.5    30.0
 3  DM_Low	           6.0    30.0
 4  Collimator	       1.2    30.0     // NOT PRESENT in spreadsheet!
 5  DM_DSM             0.0    30.0     // this need to be present for fitting? 30 arbitrary
 6  LGS_Dichroic	    -1.9    47.0  *  // was 15nm changed on Feb 23, 2026
 7  ADC1	            -3.3     6.4  *
 8  ADC2	            -4.4     6.4  *
 9  K-Mirr1           -8.3     6.6  *
 10 K-Mirr2          -12.4     6.6  *
 11 K-Mirr3          -16.5     6.6  *
 12 NGS_Dichroic	   -23.9     6.9  *  // in reflection
 13 SCI_Objective    -29.9    48.0  *
 14 SCI_Fold_Mirr_1  -36.2     5.0  *
end
// 15 SCI_Filter 1	    -425.4    8    // can't put, create numerical error cause beam is too small
// 16 SCI_Filter 2	    -479      8    // same
// 17 Cryostat window   n/a      n/a  // could be included later but likely irrelevant as very small (WFE)
*/

func read_optics_data(file,&name,&alt,&wfe)
{
  d = rdfile(file,);
  from = where(strmatch(d,"start"))(1);
  to = where(strmatch(d,"end"))(1);
  d = rdcols(file,nskip=from,nlines=to-from-1);
  name = *d(2);
  alt = *d(3);
  wfe = *d(4);
}

status = read_optics_data("optics_data_alt_wfe.i",opt_name,opt_alt,opt_wfe);
