#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._

object HelloWorld {
    println("Beginning workflow construction")
    val helloWorld = BashEvaluate(script = "echo Hello world!")
    println("Ending workflow construction")
}