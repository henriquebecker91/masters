Gem::Specification.new do | s |
  s.name        = 'batch_experiment'
  s.version     = '1.0.2'
  s.licenses    = ['Unlicense']
  s.summary     = 'A ruby script that distributes system commands between cpu cores, and save their output.'
  s.description = ''
  s.author      = 'Henrique Becker'
  s.email       = 'henriquebecker91@gmail.com'
  s.files       = [
    'lib/batch_experiment.rb',
    'lib/batch_experiment/extractor.rb',
    'lib/batch_experiment/sample_extractors.rb',
    'examples/sample_batch.rb',
    'examples/ukp_batch.rb',
    'examples/bible.txt',
    'examples/taoteching.txt',
    'examples/example_batch.rb',
    'examples/debug_batch.rb',
    'README.md',
  ]
  s.add_runtime_dependency 'childprocess', '~> 0.5'
  s.homepage    = 'https://rubygems.org/gems/batch_experiment'
end
