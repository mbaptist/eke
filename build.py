#!/usr/bin/env python

import os
import subprocess
import argparse

def Main():

	parser = build_parser()
	args = parser.parse_args()

	cmake_cmd = "/usr/bin/cmake"
	
	top_dir = os.path.dirname(os.path.realpath(__file__))
	
	if (args.build_prefix is not None):
		if (args.build_prefix=="build"):
			build_dir = os.path.join(top_dir,args.build_prefix)
		else:
			build_dir = os.path.abspath(args.build_prefix)	

	if (args.install_prefix is not None):
		if (args.install_prefix=="install"):
			install_dir = os.path.join(top_dir,args.install_prefix)
		else:     
			install_dir = os.path.abspath(args.install_prefix)

	if (args.blitz_prefix is not None):
		if (args.blitz_prefix=="install"):
			blitz_dir = os.path.join(top_dir,args.blitz_prefix)
		else:     
			blitz_dir = os.path.abspath(args.blitz_prefix)
	
	if(not os.path.exists(blitz_dir)):
	  print "Please provide the location of an existing blitz (<=0.9) installation"
	  sys.exit(1)

	print build_dir
	try:
		os.makedirs(build_dir)
	except:
		pass	
	os.chdir(build_dir)

	cmdline = [cmake_cmd]
	cmdline.append("-DCMAKE_INSTALL_PREFIX="+install_dir)
	cmdline.append("-DBLITZ_PREFIX="+blitz_dir)
	cmdline.append(top_dir)
	os_run(cmdline)

	os_run(["make","all"])

	try:
		os.makedirs(install_dir)
	except:
		pass
	os_run (["make","install"])


def os_run(cmd_line):
	p=subprocess.Popen(cmd_line)
	p.wait()

def build_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--build_prefix', dest='build_prefix', metavar='PATH', default="build",
		  help='path for building (default: ./build)')
	parser.add_argument('--install_prefix', dest='install_prefix',metavar='PATH', default="install",
		  help='installation prefix (default: ./install')
	parser.add_argument('--blitz_prefix', dest='blitz_prefix', metavar='PATH', default="install",
		  help='path of Blitz++ installation prefix (default: ./install)')
	return parser


if __name__ == "__main__":
	Main()

