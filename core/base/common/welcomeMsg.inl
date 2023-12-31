// TTK Welcome message.
// Julien Tierny <julien.tierny@sorbonne-universite.fr>
// January 2020.

// "TTK                   (c) 2024"

printMsg(
  debug::output::BOLD
    //+ " _____ _____ _  __                    __  __    ____   ___ ____  _____"
    + " _____ _____ _  __                   __  __    ____   ___ ____  _  _"
    + debug::output::ENDCOLOR,
  debug::Priority::PERFORMANCE,
  debug::LineMode::NEW,
  stream);
printMsg(debug::output::BOLD
    //+ "|_   _|_   _| |/ /                   / /__\\ \\  |___ \\ / _ " "\\___ \\|___ /"
    + "|_   _|_   _| |/ /                  / /__\\ \\  |___ \\ / _ \\___ \\| || |"
           + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);
printMsg(
  debug::output::BOLD
    //+ "  | |   | | | ' /                   | |/ __| |   __) | | | |__) | |_ \\"
    + "  | |   | | | ' /                  | |/ __| |   __) | | | |__) | || |_"
    + debug::output::ENDCOLOR,
  debug::Priority::PERFORMANCE,
  debug::LineMode::NEW,
  stream);
printMsg(
  debug::output::BOLD
//    + "  | |   | | | . \\                   | | (__| |  / __/| |_| / __/ ___) |"
    + "  | |   | | | . \\                  | | (__| |  / __/| |_| / __/|__   _|"
    + debug::output::ENDCOLOR,
  debug::Priority::PERFORMANCE,
  debug::LineMode::NEW,
  stream);
printMsg(debug::output::BOLD
//    + "  |_|   |_| |_|\\_\\                  | |\\___| | " "|_____|\\___/_____|____/"
    + "  |_|   |_| |_|\\_\\                 | |\\___| | |_____|\\___/_____|  |_|"
           + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);
printMsg(debug::output::BOLD + "                                    \\_\\  /_/"
           + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);

printMsg(debug::output::BOLD + "Welcome!" + debug::output::ENDCOLOR,
         debug::Priority::PERFORMANCE,
         debug::LineMode::NEW,
         stream);
