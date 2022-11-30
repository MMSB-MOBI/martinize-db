import Express from 'express';


export enum ErrorType {
  /** SERVER ERRORS */
  Server = 1,

  /** CLIENT ERRORS: Not Found */
  NotFound = 101,
  UserNotFound,
  ElementNotFound,

  /** CLIENT ERRORS: Unauthorized access */
  Forbidden = 201,
  TokenInvalid,
  InvalidPassword,
  Unallowed,
  UserNotApproved,

  /** CLIENT ERRORS: Bad Request */
  Format = 301,
  MissingParameters, 
  UsernameExists,
  EmailExists,
  InvalidMethod,
  TooManyFiles,
  FileTooLarge,
  InvalidMoleculeFiles,
  MissingFiles,

  /** When add molecule errors */
  UnknownParent,
  InvalidName,
  InvalidAlias,
  InvalidVersion,
  InvalidCategory,
  NameAlreadyExists,
  AliasAlreadyExists,
  VersionAlreadyExists,
  InvalidMartinizeVersion,
  InvalidForceField,
  MoleculeNotFound,

  /** When user updates his informations */
  InvalidUsername,
  InvalidEmail,

  IncorrectItpName,
  MissingTopFiles,
  ConsistencyVersionTree,
  HistoryNotFound,
  HistoryFilesNotFound, 

  /** MARTINIZE Errors */
  MartinizeRunFailed = 401,
  MartinizeNoOutput, //402
  ContactMapFailed, //403
  GOComputationFailed, //404
  ElasticNetworkFailed, //405

  /** JM Error */
  JMError = 501,
  JobNotFound, 
  UserNotProvided,
  JobNotProvided
}

const ErrorsToText = {
  [ErrorType.Server]: [500, "Server error"],
  [ErrorType.NotFound]: [404, "Page not found"],
  [ErrorType.UserNotFound]: [404, "User not found"],
  [ErrorType.ElementNotFound]: [404, "Element not found"],
  [ErrorType.MoleculeNotFound]: [404, "You can't edit a molecule that does not exists"],
  [ErrorType.Forbidden]: [403, "Access restricted"],
  [ErrorType.Unallowed]: [403, "You don't have the right to do that"],
  [ErrorType.TokenInvalid]: [403, "Invalid or expired token"],
  [ErrorType.InvalidPassword]: [403, "Invalid password"],
  [ErrorType.UserNotApproved]: [403, "Your account has not been approved yet"],
  [ErrorType.Format]: [400, "Parameter format is invalid"],
  [ErrorType.MissingParameters]: [400, "Missing parameters"],
  [ErrorType.MissingFiles]: [400, "Missing files attached to request, at least one ITP file and one PDB file is required"],
  [ErrorType.TooManyFiles]: [400, "Too many files specified"],
  [ErrorType.FileTooLarge]: [400, "Sended file is too large"],
  [ErrorType.UnknownParent]: [400, "Unknown molecule parent"],
  [ErrorType.InvalidCategory]: [400, "Invalid category"],
  [ErrorType.InvalidName]: [400, "Invalid molecule name"],
  [ErrorType.NameAlreadyExists]: [400, "Molecule name already taken"],
  [ErrorType.InvalidAlias]: [400, "Invalid molecule alias"],
  [ErrorType.AliasAlreadyExists]: [400, "Molecule alias already taken"],
  [ErrorType.InvalidVersion]: [400, "Invalid molecule version"],
  [ErrorType.VersionAlreadyExists]: [400, "Molecule version already taken"],
  [ErrorType.InvalidMoleculeFiles]: [400, "Sended molecule files (ITP, PDB/GRO) are incorrect"],
  [ErrorType.InvalidForceField]: [400, "Unknown force field"],
  [ErrorType.InvalidMartinizeVersion]: [400, "Martinize version does not exists"],
  [ErrorType.InvalidUsername]: [400, "Username contains illegal characters or is less than 2 characters long"],
  [ErrorType.InvalidEmail]: [400, "Email address is invalid"],
  [ErrorType.UsernameExists]: [409, "Username already exists"],
  [ErrorType.EmailExists]: [409, "Email already exists"],
  [ErrorType.InvalidMethod]: [405, "Method not allowed"],
  [ErrorType.MartinizeRunFailed]: [400, "Martinize run failed"],
  [ErrorType.MartinizeNoOutput] : [404, "No output created by martinize"],
  [ErrorType.ContactMapFailed] : [404, "Computation of contact map failed"],
  [ErrorType.GOComputationFailed] : [404, "Computation of GO model failed"],
  [ErrorType.ElasticNetworkFailed] : [404, "Computation of elastic network failed"],
  [ErrorType.IncorrectItpName]: [400, "The itp file name could not be parsed, please check the syntax"],
  [ErrorType.MissingTopFiles]: [400, "Missing files attached to request, there must be one top file for each itp"],
  [ErrorType.JMError] : [400, "Error with Job manager"],
  [ErrorType.HistoryNotFound] : [404, "History not found", "HistoryNotFound"], 
  [ErrorType.UserNotProvided] : [400, "User not provided to server"], 
  [ErrorType.JobNotProvided] : [400, "Job not provided to server"], 
  [ErrorType.HistoryFilesNotFound] : [404, "Job result files not found on distant file system", "HistoryFileNotFound"], 
  [ErrorType.JobNotFound] : [404, "Job doesn't exist in database", "JobNotFound"], 
  [ErrorType.ConsistencyVersionTree] : [400, "Several tree version for this molecule", "ConsistencyVersionTree"]
  
};

export default new class Errors {
  send(code: ErrorType, res: Express.Response, additionnal?: object) {
    if (res.headersSent) {
      return this.throw(code);
    }

    const [http_code, message] = ErrorsToText[code] as [number, string];

    res.status(http_code).json({
      code,
      message,
      ...(additionnal || {})
    });
  }

  /**
   * Throw the given error.
   * If you're in a promise, make sure the error is correctly catched and sent !
   */
  throw(code: ErrorType, additionnal?: object) : never {
    console.log("THROW ERROR")
    const [http_code, message] = ErrorsToText[code] as [number, string];
    console.log("http", "message", http_code, message)
    console.log("code", code)
    console.log(additionnal)
    throw new ApiError(String(http_code), {
      code,
      message,
      ...(additionnal || {})
    });
  }

  /**
   * Make the given error.
   * If you're in a promise, make sure the error is correctly catched and sent !
   */
  make(code: ErrorType, additionnal?: object) : ApiError {
    const [http_code, message, type] = ErrorsToText[code] as [number, string, string?];

    return new ApiError(String(http_code), {
      code,
      message,
      type,
      ...(additionnal || {})
    });
  }
};

export class ApiError extends Error {
  code: ErrorType;

  constructor(
    public message: string, 
    public data: { code: ErrorType, message: string, type?: string, [additionnal: string]: any }
  ) { 
    super(message); 
    this.name = "ApiError";
    this.code = this.data.code;
  }

  toJSON() {
    return {
      error: this.name,
      http_code: this.message,
      code: this.code,
      message: this.data.message,
      data: this.data,
      type : this.data.type, 
    };
  }
}

export function isCouchNotFound(e : any){
  return e.scope === "couch" && e.error === "not_found"
}

export function notFoundOnFileSystem(e: any){
  return e.code === "ENOENT"
}
