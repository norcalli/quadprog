local vim = vim
local line_splitter = require 'line_splitter'
local nvu = require 'nvim_utils'
local format = string.format
local concat = table.concat
local startswith = require 'startswith'
local insert = table.insert
local remove = table.remove


local patterns = [===[
s/const Vector<double>& \([a-zA-Z_]\+\)/\1: \&Vector
s/const Matrix<double>& \([a-zA-Z_]\+\)/\1: \&Matrix
s/Vector<double>& \([a-zA-Z_]\+\)/\1: \&mut Vector
s/Matrix<double>& \([a-zA-Z_]\+\)/\1: \&mut Matrix
s/\[\([^\]]\+\)\]\[\([^\]]\+\)\]/[(\1, \2)]/g
s/\(Vector<double>\)\s\+\(\w\+\)(\([^)]\+\))\s*[,;]\?\s*/let mut \2 = Vector::zeros(\3);\r\1 /g
s/\(Vector<int>\)\s\+\(\w\+\)(\([^)]\+\))\s*[,;]\?\s*/let mut \2 = IVector::zeros(\3);\r\1 /g
s/\(Vector<bool>\)\s\+\(\w\+\)(\([^)]\+\))\s*[,;]\?\s*/let mut \2 = BVector::default(\3);\r\1 /g
s/\(double\)\s\+\(\([a-zA-Z_0-9]\+\),\?\)\+/let mut \3: f64 = 0.0f64;\r\1 /g
s/\(unsigned int\)\s\+\(\([a-zA-Z_0-9]\+\),\?\)\+/let mut \3: usize = 0;\r\1 /g
s/for (\([ijk]\) = \([^;]\+\); \1 <\s*\([^;]\+\); \1++)/for \1 in (\2)..(\3)/
s/for (\([ijk]\) = \([^;]\+\); \1\s*>=\s*\([^;]\+\);\s*\1--)/for \1 in ((\3)..(\2 + 1)).rev()/
s/scalar_product(\([^,)]\+\),\s*\([^,)]\+\))/\1.dot(\&\2)/g
s/std::max<double>(\([^,)]\+\),\s*\([^,)]\+\))/\1.max(\2)/g
s/std::min(\([^,)]\+\),\s*\([^,)]\+\))/\1.min(\2)/g
s/sqrt(\([^,)]\+\))/(\1).sqrt()/g
s/fabs(\([^,)]\+\))/(\1).abs()/g
s/std::numeric_limits<double>::epsilon()/f64::EPSILON/g
s/\([a-zA-Z0-9_]\+\)++;/\1 += 1;/
s/\([a-zA-Z0-9_]\+\)--;/\1 -= 1;/
s/goto\s\+\([a-zA-Z0-9_]\+\);/continue '\1;/
]===]

return {
  on_attach = function(fname, root)
    if vim.bo.ft == 'rust' then
      nvu.nvim_create_commands({
        Transform1 = {
          function()
            nvu.nvim_buf_transform_region_lines(0, '<', '>', VISUAL_MODE.line, function(lines)
              local r = {}
              local deleting = false
              local gotos = {}
              for _, line in ipairs(lines) do
                if not deleting and line:find("#ifdef", 1, true) then
                  deleting = true
                end
                if not deleting then
                  line = line:gsub("^%s*([a-z0-9A-Z_]+):", function(label)
                    if label == "std" then
                      return
                    end
                    gotos[#gotos+1] = label
                    r[#r+1] = format("'%s: loop {", label)
                    return ""
                  end)
                  r[#r+1] = line
                end
                if deleting and line:find("#endif", 1, true) then
                  deleting = false
                end
              end
              while #gotos > 0 do
                local label = remove(gotos)
                r[#r+1] = format("break '%s;", label)
                r[#r+1] = "}"
              end
              return r
            end)
          end;
          range = true;
          buffer = true;
        };
        Transform = {
          function()
            for pattern in line_splitter().drain(patterns) do
              for i = 1, 10 do
                vim.cmd("silent! '<,'>".."+"..i..pattern)
              end
            end
          end;
          range = true;
          buffer = true;
        };
      })
    end
  end
}
