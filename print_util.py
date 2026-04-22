import sys
import io
sys.stdout = io.TextIOWrapper(
    sys.stdout.buffer, 
    encoding="utf-8", 
    line_buffering=True
)

SUB_CHAR_IN = "aehijklmnoprstuvx0123456789"
SUB_CHAR_LIST = "ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ₀₁₂₃₄₅₆₇₈₉"

SUP_CHAR_IN =   "abcdefghijklmnoprstuvwxyzABDEGHIJKLMNOPRTUW0123456789"
SUP_CHAR_LIST = "ᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᴬᴮᴰᴱᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾᴿᵀᵁᵂ⁰¹²³⁴⁵⁶⁷⁸⁹"

SUB_MAP = str.maketrans(SUB_CHAR_IN, SUB_CHAR_LIST)
SUP_MAP = str.maketrans(SUP_CHAR_IN, SUP_CHAR_LIST)

def sub(n) : return str(n).translate(SUB_MAP)
def sup(n) : return str(n).translate(SUP_MAP)
