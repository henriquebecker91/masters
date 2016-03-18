Gem::Specification.new do | s |
  s.name        = 'batch_experiment'
  s.version     = '0.1.0'
  s.licenses    = ['Public Domain', 'Unlicense']
  s.summary     = 'A ruby script that distributes system commands between cpu cores, and save their output.'
  s.description = ""
  s.author      = 'Henrique Becker'
  s.email       = 'henriquebecker91@gmail.com'
  s.files       = [
    'lib/batch_experiment.rb',
    'lib/batch_experiment/extractor.rb',
    'lib/batch_experiment/sample_extractors.rb',
    'examples/sample_batch.rb',
    'examples/ukp_batch.rb',
  ]
  s.add_runtime_dependency 'childprocess', '~> 0.5'
  s.homepage    = 'https://rubygems.org/gems/batch_experiment'
end
