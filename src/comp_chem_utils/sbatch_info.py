#!/usr/bin/env python
__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"


class sbatch_option(object):
   """This class is used to deal with sbatch options"""

   def __init__(self):
      self.command = ''
      self.value = ''
      self.explanation = ''

   def output(self):
      return '{}={}'.format(self.command, self.value)

   def detailed_output(self):
      return [
            'Option: {}'.format(self.command),
            'Current value: {}'.format(self.value),
            'Explanation: {}'.format(self.explanation)
      ]

   def modify_from_user(self):
      print('')
      for line in self.detailed_output():
         print(line)
      nv = raw_input('New value (empty cancels):')
      if nv:
         self.value = nv

   def set_default(self, key):
      self.command = '--'+key

      # ADD NEW OPTIONS HERE:
      if key=='time':
         self.set( '0-1:00' )
         self.explanation = 'Max wall time, format: days-hours:minutes'

      elif key=='mem':
         self.set( '120G' )
         self.explanation = 'e.g. 1G (default unit is MB can be change with suffix K,M,G,T)'

      elif key=='partition':
         self.set( 'debug' )
         self.explanation = 'debug or parallel'

      elif key=='nodes':
         self.set( '1' )
         self.explanation = 'Number of nodes to use'

      #elif key=='ntasks':
      #   # should not be set directly, instead set nodes and cpus-per-task
      #   self.set( '28' )
      #   self.explanation = 'Number of MPI processes (default is 28)'

      elif key=='ntasks-per-node':
         self.set( '28' )
         self.explanation = 'Number of MPI processes per node (default is 28, 1 per core)'

      elif key=='job-name':
         self.set( '' )
         self.explanation = 'Job name'


   def set(self, value):
      self.value = value


# ADD NEW OPTIONS HERE:
options = [
   'time',
   'mem',
   'partition',
   'nodes',
   'ntasks-per-node',
   'job-name'
]

sbatch_info = {}

for option in options:
   sbatch_info[option] = sbatch_option()

for key, value in sbatch_info.items():
   value.set_default(key)


def sbatch_line(script_name):
   """concatenate sbatch info into a submission line"""

   line = ['sbatch']

   for key in sbatch_info:
      line.append( sbatch_info[key].output() )

   line.append( script_name )

   return ' '.join( line )

