module Extractor
  # For when there's a field whose value is after '<field>: '.
  def self.get_field(lines, field)
    lines.grep(/^#{field}: .*/).each { | l | return l.match(/:[\t ]+(.*)/)[1] }
    ''
  end

  # For when there's a field whose value is in the next line.
  def self.get_hfield(lines, field)
    if ix = lines.find_index(field) then lines[ix + 1] else '' end
  end

  # Return the field names for each of the elements returned by
  # extract. Ex.: ['Time', 'Max Mem Use', 'opt', ... ]
  def names
    fail 'This method should have been overwritten by a subclass.'
  end

  def extract(content)
    extract_from_lines(content.lines.map! { | l | l.chomp! })
  end

  # Extract an array of values from the command output. This array has the same
  # size as the one returned by field_names.
  def extract_from_lines(lines)
    fail 'This method should have been overwritten by a subclass.'
  end
end

