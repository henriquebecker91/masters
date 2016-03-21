# Module that defines the interface used for extracting info from other
# programs output. You don't need to include it in your object, will suffice
# that the object (that you will use to extract info from the output) has the
# ::names and ::extract methods defined.
module Extractor
  # Find a line in the following format: "field: value", return value. 
  #
  # @param lines [Array<String>] Program output, broken in lines.
  # @param field [String] String to be found at the lines in the following
  #   pattern: 'field: value'.
  #
  # @return [String] The 'value' as a string or, if 'field' isn't found, an
  #   empty string.
  def self.get_field(lines, field)
    lines.grep(/^#{field}: .*/).each { | l | return l.match(/:[\t ]+(.*)/)[1] }
    ''
  end

  # Find a line in the following format: "field<linebreak>value", return value.
  #
  # @param lines [Array<String>] Program output, broken in lines.
  # @param field [String] String to be found at the lines in the following
  #   pattern: 'field<lineabreak>value'. Only include the linebreak at the
  #   field, if the lines array have linebreaks at the end of every string
  #   element.
  #
  # @return [String] The 'value' as a string or, if 'field' isn't found, an
  #   empty string.
  def self.get_hfield(lines, field)
    if ix = lines.find_index(field) then lines[ix + 1] else '' end
  end

  # Return the field names for each of the elements returned by extract. Ex.:
  #   ['Time', 'Max Mem Use', 'opt', ... ]
  # 
  # @note To be on the safe side you should create a new array at each call. If
  #   you always return a reference to the same array the array can be
  #   modified.
  #
  # @return [Array<String>] The strings that will be used to make the column
  #   names at the BatchExperiment.experiment method.
  def names
    fail 'This method should have been overwritten by a subclass.'
  end

  # Extract N values of some program output, where N is equal to #names.size.
  # If some value isn't present, the array should remain the same size, and the
  # value at the corresponding position should be nil.
  # 
  # @param content [String] Program output.
  # @return [Array<String>] The N extracted values, as strings.
  def extract(content)
    extract_from_lines(content.lines.map! { | l | l.chomp! })
  end

  # Optionally, you can define this method instead of #extract. The #extract
  # method will call this method if not overrided.
  #
  # @param lines [Array<String>] Program output, broken in lines,
  #   and the line string elements don't end in linebreak.
  # @return [Array<String>] The N extracted values, as strings.
  def extract_from_lines(lines)
    fail 'This method should have been overwritten by a subclass.'
  end
end

