#! /usr/bin/env python 

'''
This is a small set of routines that creates html files


Command line usage (if any):

	usage: xhtml.py  will run a basic test of the routines

Description:  

	Each routine returns a string.  As one builds up an html page
	one normally, adds the string that is returnted to a pre-existing
	string.  When finished one writes the accumulated string to a file

Primary routines:

Notes:
									   
History:

110505 ksl Coding begun
161031 ksl Renamed to xhtml.py to avoid conflict with anaconda

'''

import sys

def begin(title='Title'):
	'''
	Start an html page 
	'''

	string='''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html lang="en">
<head>
	<meta http-equiv="content-type"
	content="text/html; charset=ISO-8859-1">
<title>%s</title>
</head>
	<body>
<h1>%s</h1>
''' % (title,title)
	return string

def end():
	'''
	Finish and html page
	'''

	string='''\n</body> \n</html>'''
	return string

def paragraph(content='This is a paragraph.  Who knows whether to continue'):
	'''
	Add a paragraph to the page
	'''
	string = '\n<p> \n	%s \n</p>\n' % content

	return string

def link(text='here',href='foo.txt'):
	'''
	Make a link to the page

	Note: generally speaking one needs to put this in a string for a paragraph
	'''
	string='''
	<a
	 href="%s">%s
	 </a>''' % (href,text)

	return string

def image(image='test.png',alt='Thumbnails',width=400,height=400):
	'''
	Add a centered image to a page.  If height is 0, it is left out, so the
	image keeps its own aspect ratio at the given width.
	'''
	string='''\n<div style="text-align: center;">'''
	if height>0:
		string=string+'''\n	<img alt="%s" src="%s" width="%d" height="%d"> \n''' % (alt,image,width,height)
	else:
		string=string+'''\n	<img alt="%s" src="%s" width="%d"> \n''' % (alt,image,width)
	string=string+'''</div>\n'''

	return string 

def table(lines):
	'''
	Format a table in html.  Lines should be a set of records, with each record containing the
	values for one row.  
	'''

	string='''\n<table border="1" cellpadding="2" cellspacing="2" width="100%">\n'''
	for line in lines:
		row='<tr>\n'
		for word in line:
			row=row+'	<td> %s </td>\n' % word
		row=row+'<tr>\n'
		string=string+row
	trailer='''
</tbody>
</table>
'''
	string=string+trailer 
	return string

def hline(size=3,width=100):
	'''
	Draw a horizontal line with thickness given by size and horizontal width given
	by width in percent
	'''

	string='''\n<hr size="%d" width="%d'''  % (int(size),int(width))
	string=string+'%">\n' # Note that this funny construction is needed to handle the % symbol

	return string

def h1(line):
	'''
	Add a header line
	'''
	string='''<h1>%s</h1>''' % line
	return string

def h2(line):
	'''
	Add a header line
	'''
	string='''<h2>%s</h2>''' % line 
	return string

def h3(line):
	'''
	Add a header line
	'''
	string='''<h3>%s</h3>''' % line
	return string

def add_list(lines):
	'''
	Add a simple list

	Note - This routine is named add_list rather than list
	to avoid confusion with python lists
	'''
	string='\n<ul>'
	for line in lines:
		line=line.strip()
		string=string+'\n	<li>%s</li>' % line
	string=string+'\n</ul>'

	return string

def preformat(lines):
	'''
	Add a set of preformated lines, that is text.
	'''

	string='<pre>'
	for line in lines:
		line=line.replace('\n','')
		string=string+line+'<br>'
	string=string+'</pre>'
	return string

def test(filename='test.html'):
	'''
	Test the various html subroutines
	'''

	string=begin('Hi Knox')
	string=string+paragraph('Tell it all, Knox')

	link_string=link('here','ds9.py')

	mypar='This file can be found %s' % link_string

	string=string+paragraph(mypar)



	line=[['test1',4],['test2',4]]
	string=string+table(line)

	string=string+hline()

	string=string+image('file:///data/Projects/Persistence/Test/D/Visit19/Persist/Figs/ib6v19cpq_subtract.png')

	string=string+hline()

	line=['Hi Knox','Goodbye Knox']

	string=string+add_list(line)

	lines=['this is test of preformating\n','hello','world','tomorrow is another day']

	string=string+preformat(lines)

	string=string+end()


	g=open(filename,'w')
	g.write(string)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		print('usage: xhtml.py test.html')
