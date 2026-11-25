#!/usr/bin/env python3
"""Reproduce and verify the complete independent check of the 20 Niemeier types."""
import ast, pathlib, re, subprocess, sys, os
ROOT=pathlib.Path(__file__).resolve().parent

def run(script,log,*args):
 path=ROOT/log
 with path.open('w') as f:
  subprocess.run([sys.executable,'-u',str(ROOT/script),*args],stdout=f,stderr=subprocess.STDOUT,check=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
 return path.read_text().splitlines()

def first(lines):
 records={};data={};i=None
 for l in lines:
  if l.startswith('START '):i=int(l.split()[1]);assert i not in records
  if l.startswith('RESULT '):
   m=re.fullmatch(r'RESULT (\d+) max bound (\d+) explicit markings (\d+) remaining (\d+) passed bound (\d+)',l);assert m
   k,b,c,n,a=map(int,m.groups());assert k==i and a<=144;data[i]=(n,a)
  if l.startswith('REMAINING RECORDS '):
   records[i]=ast.literal_eval(l.split(' ',2)[2]);assert len(records[i])==data[i][0]
 assert set(records)==set(data)==set(range(20))
 return records,data

def second(lines,records):
 ans={}
 for l in lines:
  if l.startswith('RESULT '):
   m=re.fullmatch(r'RESULT (\d+) records (\d+) max axes bound (\d+) remaining (\d+) max base (\d+) deletions (\d+)',l);assert m
   i,n,b,remain,base,d=map(int,m.groups());assert i not in ans and n==len(records[i]) and remain==0 and b<=72
   ans[i]=(b,base,d)
 assert set(ans)=={i for i,r in records.items() if r}
 return ans

if __name__=='__main__':
 a=run('niemeier20_independent.py','niemeier20_first_pass.txt');records,data=first(a)
 print('PASS: first pass, 20 types,',sum(map(len,records.values())),'marked exceptions',flush=True)
 b=run('niemeier20_global.py','niemeier20_second_pass.txt',str(ROOT/'niemeier20_first_pass.txt'));ans=second(b,records)
 print('PASS: every marked exception excluded',flush=True)
 theta=run('niemeier20_theta_check.py','niemeier20_shell_counts.txt');assert [int(l.split()[0]) for l in theta]==list(range(20))
 print('PASS: all norm 2,4,6 shell counts checked',flush=True)
 for i in range(20):
  print(i,'exceptions',data[i][0],'final signed-point bound',max(data[i][1],2*ans.get(i,(0,0,0))[0]),'deletions',ans.get(i,(0,0,0))[2],flush=True)
 print('ALL 20 TYPES VERIFIED',flush=True)
