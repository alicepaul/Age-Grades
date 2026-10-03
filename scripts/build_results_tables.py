import csv
from pathlib import Path
root=Path(__file__).resolve().parents[1]/'Output'
dest=root/'latex'
dest.mkdir(exist_ok=True)
def read(name):
 with (root/(name+'.csv')).open() as f: return list(csv.DictReader(f))
def num(v): return f'{float(v):.2f}'
def pv(v):
 x=float(v)
 if x<.001:
  return r'$<0.001$'
 return f'{x:.4f}'
def ci(r): return '$['+num(r['ci_lower'])+', '+num(r['ci_upper'])+']$'
def table(name,source,caption,headers,rows,note):
 s='% Source: Output/'+source+'.csv\n% Requires booktabs.\n'
 s+='\\begin{table}[htbp]\n\\centering\n\\small\n'
 s+='\\caption{'+caption+'}\n\\label{tab:'+name.replace('_','-')+'}\n'
 s+='\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{'+('l'*2+'r'*(len(headers)-2))+'}\n\\toprule\n'
 s+=' & '.join(headers)+r' \\'+'\n\\midrule\n'
 for row in rows: s+=' & '.join(row)+r' \\'+'\n'
 s+='\\bottomrule\n\\end{tabular}\n\\par\\medskip\n\\begin{minipage}{\\textwidth}\n\\footnotesize\n'+note+'\n\\end{minipage}\n\\end{table}\n'
 (dest/(name+'.tex')).write_text(s);return s
ranknote=r'Gaussian identity-link GEEs use working independence and race-clustered robust sandwich standard errors. Tests are two-sided Wald tests; intervals are robust Wald 95\% confidence intervals for the same mean contrasts. Matched winners use an intercept-only model for female-minus-male differences within race--age-group pairs; all finishers use a sex indicator with male as reference. $n$ counts pairs or runners, respectively; clusters are 216 races. Each pair or runner has equal weight. Differences are percentage points. Inference assumes independent races and conditions on the fitted standards.'
mar=read('marathon_test_table'); met=read('method_comparison')
for r in met:
 for standard, column in [('Existing', 'existing_gap'),
                          ('Proposed', 'proposed_gap')]:
  other=next(x for x in mar if x['Standard']==standard and
             x['Analysis']==r['Analysis'])
  assert abs(float(r[column])-float(other['mean_difference'])) < 1e-9
 assert abs(float(r['change_in_mean_gap']) -
            (float(r['proposed_gap'])-float(r['existing_gap']))) < 1e-9
parts=[]
parts.append(table('method_comparison','method_comparison','Direct comparison of existing and proposed age-grading standards.', ['Analysis','$n$','Existing gap','Proposed gap','Change','95\\% CI','GEE $p$'], [[r['Analysis'],f"{int(r['n']):,}",num(r['existing_gap']),num(r['proposed_gap']),num(r['change_in_mean_gap']),ci(r),pv(r['p_value'])] for r in met],ranknote+r' Change is the proposed gap minus the existing gap. Direct-comparison models use within-runner changes in age grade, paired across sex for matched winners.'))
parts.append(table('marathon_test_table','marathon_test_table','Sex differences in marathon age grades.', ['Standard','Analysis','$n$','Mean gap','95\\% CI','GEE $p$'], [[r['Standard'],r['Analysis'],f"{int(r['n']):,}",num(r['mean_difference']),ci(r),pv(r['p_value'])] for r in mar],ranknote))
# The requested best_test_table.csv is absent; records contain single-age bests.
records=read('records_test_table')
parts.append(table('best_test_table','records_test_table','Sex differences in age grades for single-age best performances.', ['Standard','Event','Pairs','Female mean','Male mean','Mean gap','95\\% CI','$p$'], [[r['Standard'],r['Event'],r['n_pairs'],num(r['mean_female']),num(r['mean_male']),num(r['mean_difference']),ci(r),pv(r['p_value'])] for r in records],r'Means are age-grade percentages; gaps are female minus male in percentage points. Tests are two-sided paired $t$-tests, matching by age within event. Intervals are nominal 95\% paired $t$ intervals. Source: \texttt{records\_test\_table.csv}, used for single-age bests because \texttt{best\_test\_table.csv} was not present.'))
preview=r'''\documentclass[11pt]{article}
\usepackage[margin=0.65in]{geometry}
\usepackage{booktabs}
\begin{document}
'''+ '\n\\clearpage\n'.join(parts)+'\n\\end{document}\n'
(dest/'results_tables.tex').write_text(preview)
