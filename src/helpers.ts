import logger from "./logger";

export function isDebugMode() {
  return logger.level === "debug" || logger.level === "silly";
}
