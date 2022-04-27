// TTK Welcome message.
// Julien Tierny <julien.tierny@sorbonne-universite.fr>
// January 2020.

// "TTK     (c) 2022"

printMsg(
  debug::output::BOLD
    + " _____ _____ _  __                    __  __    ____   ___ ____  ____"
    + debug::output::ENDCOLOR,
  debug::Priority::PERFORMANCE,
  debug::LineMode::NEW,
  stream);
printMsg(debug::output::BOLD
           + "|_   _|_   _| |/ /                   / /__\\ \\  |___ \\ / _ "
             "\\___ \\|___ \\"
           + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);
printMsg(
  debug::output::BOLD
    + "  | |   | | | ' /                   | |/ __| |   __) | | | |__) | __) |"
    + debug::output::ENDCOLOR,
  debug::Priority::PERFORMANCE,
  debug::LineMode::NEW,
  stream);
printMsg(
  debug::output::BOLD
    + "  | |   | | | . \\                   | | (__| |  / __/| |_| / __/ / __/"
    + debug::output::ENDCOLOR,
  debug::Priority::PERFORMANCE,
  debug::LineMode::NEW,
  stream);
printMsg(debug::output::BOLD
           + "  |_|   |_| |_|\\_\\                  | |\\___| | "
             "|_____|\\___/_____|_____|"
           + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);
printMsg(debug::output::BOLD + "                                     \\_\\  /_/"
           + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);

printMsg(debug::output::BOLD + "Welcome!" + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);
